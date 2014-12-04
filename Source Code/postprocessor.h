#ifndef _POSTPROCESS_H_
#define _POSTPROCESS_H_

#include <fstream>
#include <queue>
#include <boost/regex.hpp>
#include <iomanip>
#include <iterator>
#include "trivial_match_map.h"
#include "pipeline.h"

typedef Candidate payload;

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


class RemoveTrivialMatchesFilter {
private:
    TrivialMatchMap trivial_matches_map;
    std::queue<Candidate> trivial_match_prevention_queue;
public:
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
                std::cout << "loc: " << std::setw(10) << oldest_candidate.query_loc << ", candidate loc: " << std::setw(10) << oldest_candidate.loc << ", dist: " << std::setw(10) << oldest_candidate.dist << std::endl;
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
    void init(std::vector<std::string> const& input_file_names, std::function<void(payload const&)> yield) {
        for (std::string const& input_file_name: input_file_names) {
            input_pq.push(new InputFile(input_file_name));
        }
        while ( ! input_pq.empty()) {
            InputFile& input_file = *(input_pq.top());
            input_pq.pop();
            size_t curr_query_loc = input_file->query_loc;
            size_t prev_query_loc;
            do {
                Candidate new_candidate = *input_file;
                yield(new_candidate);
                ++input_file;
                prev_query_loc = new_candidate.query_loc;
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
    Pipeline<Candidate, std::vector<std::string>, InputMergeProducer, RemoveTrivialMatchesFilter> pipeline;
    std::vector<std::string> input_file_names;
public:
    PostProcessor(std::vector<std::string> const& input_file_names): input_file_names(input_file_names) {}
    void run() {
        pipeline.go(input_file_names);
    }
};

#endif
