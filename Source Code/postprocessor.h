#ifndef _POSTPROCESS_H_
#define _POSTPROCESS_H_

#include <fstream>
#include <queue>
#include <boost/regex.hpp>
#include <iomanip>
#include "trivial_match_map.h"

class PostProcessor {
    class InputFile;
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

    class InputFile {
    private:
        std::ifstream input;
        bool eof;
        static boost::regex query_cand_dist_expr; //matches input file line fields
        size_t curr_query_pos;
        size_t curr_candidate_pos;
        double curr_distance;
    public:
        bool operator>(InputFile const& rhs) const {
            //compare input files by query position so they can be merged in order
            return curr_query_pos > rhs.curr_query_pos;
        }
        void next() {
            //read from the input file until we find a valid line or we run out of lines
            while ( ! eof) {
                std::string line;
                if (input >> line) {
                    std::string::const_iterator start = line.begin();
                    std::string::const_iterator end = line.end();
                    boost::smatch pieces;
                    if (regex_search(start, end, pieces, query_cand_dist_expr)) {
                        curr_query_pos = stoi(pieces[1]);
                        curr_candidate_pos = stoi(pieces[2]);
                        curr_distance = stod(pieces[3]);
                        break;
                    }
                } else
                    eof = true; //end of input file
            }
        }
        InputFile(std::string const& input_file_name): input(input_file_name), eof(false) {
            next();
        }
        size_t get_curr_query_pos() const { return curr_query_pos; }
        size_t get_curr_candidate_pos() const { return curr_candidate_pos; }
        double get_curr_distance() const { return curr_distance; }
        bool is_eof() const { return eof; };
    };

    void process_input_file_candidate(Candidate const& input_candidate, std::queue<Candidate>* trivial_match_prevention_queue, TrivialMatchMap* trivial_matches_map) {
        Candidate old_candidate;
        if (trivial_matches_map->detect_trivial_match(input_candidate, &old_candidate)
            && input_candidate.dist > old_candidate.dist)
            //found a trivial match and old candidate is better than new, so dont do anything
            return;
        //no trivial match or new candidate is better
        trivial_matches_map->push_candidate(input_candidate);
        trivial_match_prevention_queue->push(input_candidate);
        while (trivial_match_prevention_queue->size() > 0)
        {
            Candidate oldest_candidate = trivial_match_prevention_queue->front();
            if (oldest_candidate.is_query_close_to(input_candidate))
                break;
            trivial_match_prevention_queue->pop();
            Candidate candidate;
            if (trivial_matches_map->detect_trivial_match(oldest_candidate, &candidate)
                && oldest_candidate == candidate)
            {
                //candidate made it through thr trivial match checking process without being deleted
                std::cout << "loc: " << std::setw(10) << oldest_candidate.query_loc << ", candidate loc: " << std::setw(10) << oldest_candidate.loc << ", dist: " << std::setw(10) << oldest_candidate.dist << std::endl;
            }
        }
    }

public:
    PostProcessor(std::vector<std::string> const& input_file_names) {
        for (std::string const& input_file_name: input_file_names) {
            input_pq.push(new InputFile(input_file_name));
        }
    }

    void run() {
        TrivialMatchMap trivial_match_map;
        std::queue<Candidate> trivial_match_prevention_queue;
        while ( ! input_pq.empty()) {
            InputFile* input_file = input_pq.top();
            input_pq.pop();
            size_t curr_query_pos = input_file->get_curr_query_pos();
            size_t prev_query_pos;
            do {
                Candidate new_candidate;
                new_candidate.dist = input_file->get_curr_distance();
                new_candidate.query_loc = curr_query_pos;
                new_candidate.loc = input_file->get_curr_candidate_pos();
                process_input_file_candidate(new_candidate, &trivial_match_prevention_queue, &trivial_match_map);
                input_file->next();
                prev_query_pos = curr_query_pos;
                curr_query_pos = input_file->get_curr_query_pos();
            } while (prev_query_pos == curr_query_pos);
            if ( ! input_file->is_eof())
                input_pq.push(input_file);
        }
    }
};

#endif
