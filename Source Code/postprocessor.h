#ifndef _POSTPROCESS_H_
#define _POSTPROCESS_H_

#include <fstream>
#include <queue>
#include <deque>
#include <boost/regex.hpp>
#include <iomanip>
#include <iterator>
#include <map>
#include <boost/any.hpp>
#include "trivial_match_map.h"
#include "pipeline.h"
#include "subsequence.h"
#include "ucr_dtw.h"

typedef Candidate payload;

namespace PostProcess {
    class InputFile: public std::iterator<std::forward_iterator_tag, const Candidate> {
    private:
        std::ifstream input;
        bool eof;
        boost::regex query_cand_dist_expr;
        static constexpr const char* query_cand_dist_expr_str = R"(^(\d+)[ ,](\d+)[ ,](\d*\.?\d*))"; //matches input file line fields
        Candidate curr_candidate;
        void find_next_valid_line() {
            //read from the input file until we find a valid line or we run out of lines
            while ( ! eof) {
                std::string line;
                if (input >> line) {
                    std::string::const_iterator start = line.begin();
                    std::string::const_iterator end = line.end();
                    boost::smatch pieces;
                    if (regex_search(start, end, pieces, query_cand_dist_expr)) {
                        size_t curr_query_loc = stoi(pieces[1]);
                        size_t curr_candidate_loc = stoi(pieces[2]);
                        double curr_distance = stod(pieces[3]);
                        curr_candidate.query_loc = curr_query_loc;
                        curr_candidate.loc = curr_candidate_loc;
                        curr_candidate.dist = curr_distance;
                        break;
                    }
                } else
                    eof = true; //end of input file
            }
        }
    public:
        bool operator>(InputFile const& rhs) const {
            //compare input files by query position so they can be merged in order
            return curr_candidate.query_loc > rhs.curr_candidate.query_loc;
        }
        InputFile(std::string const& input_file_name): input(input_file_name), eof(false), query_cand_dist_expr(query_cand_dist_expr_str) {
            find_next_valid_line();
        }
        Candidate const& operator*() const { return curr_candidate; }
        Candidate const* operator->() const { return & operator*(); }
        Candidate const& operator++() { find_next_valid_line(); }
        bool is_eof() const { return eof; };
    };

    class CandidateOutputConsumer {
    public:
        CandidateOutputConsumer(std::map<std::string, boost::any> init_args) {}
        void go(payload const& input_candidate, std::function<void(payload const&)> yield) {
            std::cout << input_candidate.query_loc << "\t" << input_candidate.loc << "\t" << input_candidate.dist << std::endl;
        }
        void end(std::function<void()> yield) {
            yield();
        }
    };

    class SelfMatchBridgingFilter {
    private:
        std::deque<size_t> potential_self_match_query_locs;
        SubsequenceLookup& subsequences;
        UCR_DTW ucr_suite;
        size_t query_len;
        /*
            dummy_vector is passed to to the ucr_suite's dtw dist function as the lower bound subsequence (all zeros
            means no lower bound, ie don't abandon early)
        */
        std::vector<double> dummy_lb_vector;
        double get_dtw_dist(size_t query_loc1, size_t query_loc2) {
            //returns the dtw distance between the subsequences stored at the query positions pointed to by the two candidates
            Subsequence const& query1 = subsequences[query_loc1];
            Subsequence const& query2 = subsequences[query_loc2];
            double dist = ucr_suite.dtw(query1.series_normalized, query2.series_normalized, dummy_lb_vector.data(), query_len);
            return dist;
        }
    public:
        SelfMatchBridgingFilter(std::map<std::string, boost::any>& init_args):
            subsequences(*(boost::any_cast<SubsequenceLookup*>(init_args["subsequence-lookup"]))),
            ucr_suite(subsequences),
            query_len(boost::any_cast<size_t>(init_args["query-len"])),
            dummy_lb_vector(query_len, 0) {}
        void go(payload const& input_candidate, std::function<void(payload const&)> yield) {
            size_t last_added_query_loc = potential_self_match_query_locs.back();
            if (last_added_query_loc != input_candidate.query_loc) {
                //haven't already seen this query loc
                while (potential_self_match_query_locs.size() > 0)
                {
                    size_t oldest_query_loc = potential_self_match_query_locs.front();
                    Candidate oldest_query_loc_placeholder_candidate;
                    oldest_query_loc_placeholder_candidate.query_loc = oldest_query_loc;
                    if (oldest_query_loc_placeholder_candidate.is_query_close_to(input_candidate)) //xxx check this, should it be input candidate?
                        break;
                    potential_self_match_query_locs.pop_front();
                    for (size_t bridge_target_query_loc: potential_self_match_query_locs) {
                        Candidate bridge_query_loc_placeholder_candidate;
                        bridge_query_loc_placeholder_candidate.query_loc = bridge_target_query_loc;
                        if ( ! bridge_query_loc_placeholder_candidate.is_query_close_to(oldest_query_loc_placeholder_candidate))
                            //not a self-match
                            break;
                        //create a bridge between two candidates whose query positions are close
                        double dist = get_dtw_dist(oldest_query_loc, bridge_target_query_loc);
                        Candidate bridge;
                        bridge.query_loc = oldest_query_loc;
                        bridge.loc = bridge_target_query_loc;
                        bridge.dist = dist;
                        yield(bridge);
                    }
                }
                potential_self_match_query_locs.push_back(input_candidate.query_loc);
            }
            //forward the non-trivial match
            yield(input_candidate);
        }
        void end(std::function<void()> yield) {
            yield();
        }
    };

    class RemoveTrivialMatchesFilter {
    private:
        TrivialMatchMap trivial_matches_map;
        std::queue<Candidate> trivial_match_prevention_queue;
    public:
        RemoveTrivialMatchesFilter(std::map<std::string, boost::any>& init_args) {}
        void go(payload const& input_candidate, std::function<void(payload const&)> yield) {
            Candidate old_candidate;
            if (trivial_matches_map.detect_trivial_match(input_candidate, &old_candidate)
                && input_candidate.dist > old_candidate.dist)
                //found a trivial match and old candidate is better than new, so dont do anything
                return;
            //no trivial match or new candidate is better
            trivial_matches_map.push_candidate(input_candidate);
            trivial_match_prevention_queue.push(input_candidate);
            while (trivial_match_prevention_queue.size() > 0)
            {
                Candidate oldest_candidate = trivial_match_prevention_queue.front();
                if (oldest_candidate.is_query_close_to(input_candidate))
                    break;
                trivial_match_prevention_queue.pop();
                Candidate candidate;
                if (trivial_matches_map.detect_trivial_match(oldest_candidate, &candidate)
                    && oldest_candidate == candidate)
                {
                    //candidate made it through the trivial match checking process without being deleted
                    yield(oldest_candidate);
                }
            }
        }
        void end(std::function<void()> yield) {
            yield();
        }
    };

    class InputMergeProducer {
    /*
        since the motif finder generated multiple files, each containing candidates sorted by query loc, we use a priority
        queue to externally sort them, yielding candidates in order by query loc
    */
    private:
        template<typename T>
        struct deref_greater: public std::greater<T> {
            //allows us to store pointers in the priority queue and avoid excess copying
            bool operator()(const T* x, const T* y) const {
                return *x > *y;
            }
        };
        template<typename T>
        struct PointerPriorityQueue: public std::priority_queue<T*, std::vector<T*>, deref_greater<T>> {};
        PointerPriorityQueue<InputFile> input_pq;
    public:
        InputMergeProducer(std::map<std::string, boost::any>& init_args) {
            std::vector<std::string> input_file_names = boost::any_cast<std::vector<std::string>>(init_args["input-files"]);
            for (std::string const& input_file_name: input_file_names) {
                input_pq.push(new InputFile(input_file_name));
            }
        }
        void init(std::function<void(payload const&)> yield) {
            while ( ! input_pq.empty()) {
                InputFile& input_file = *(input_pq.top());
                input_pq.pop();
                size_t curr_query_loc = input_file->query_loc;
                size_t prev_query_loc;
                do {
                    Candidate const& new_candidate = *input_file;
                    yield(new_candidate);
                    prev_query_loc = new_candidate.query_loc;
                    //be sure to do everything we need to do with the candidate before incrementing the input file iterator
                    //because the const ref will change
                    ++input_file;
                    curr_query_loc = input_file->query_loc;
                } while (prev_query_loc == curr_query_loc);
                if ( ! input_file.is_eof())
                    input_pq.push(&input_file);
            }
        }
        void end(std::function<void()> yield) {
            yield();
        }
    };

    class PostProcessor {
        std::map<std::string, boost::any> init_args;
        Pipeline<Candidate,
                 InputMergeProducer,
                 RemoveTrivialMatchesFilter,
                 SelfMatchBridgingFilter,
                 CandidateOutputConsumer> pipeline;
    public:
        PostProcessor(std::vector<std::string> const& input_file_names, SubsequenceLookup* subsequences, size_t query_len):
            init_args({
                std::make_pair("input-files", input_file_names),
                std::make_pair("subsequence-lookup", subsequences),
                std::make_pair("query-len", query_len)
            }),
            pipeline(init_args) {}
        void run() {
            pipeline.go();
        }
    };
}

#endif
