#ifndef _POSTPROCESS_H_
#define _POSTPROCESS_H_

#include <fstream>
#include <queue>
#include <boost/regex.hpp>
#include <iomanip>
#include "trivial_match_map.h"

class PostProcessor {
    template<typename T>
    struct deref_greater: public std::greater<T> {
        //allows us to store pointers in the priority queue and avoid excess copying
        bool operator()(const T* x, const T* y) const {
            return *x > *y;
        }
    };
    template<typename T>
    struct PointerPriorityQueue: public std::priority_queue<T*, std::vector<T*>, deref_greater<T>> {};

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

    PointerPriorityQueue<InputFile> input_pq;
    size_t num_trivial_matches;
    size_t num_all_matches;

public:
    PostProcessor(std::vector<std::string> const& input_file_names) {
        for (std::string const& input_file_name: input_file_names) {
            input_pq.push(new InputFile(input_file_name));
        }
    }

    void process_input_file_line(InputFile const& input_file, TrivialMatchMap* trivial_matches_map) {
        size_t curr_query_pos = input_file.get_curr_query_pos();
        size_t curr_candidate_pos = input_file.get_curr_candidate_pos();
        double curr_distance = input_file.get_curr_distance();
        Candidate new_candidate;
        new_candidate.dist = curr_distance;
        new_candidate.query_loc = curr_query_pos;
        new_candidate.loc = curr_candidate_pos;
        Candidate old_candidate;
        if (trivial_matches_map->detect_trivial_match(new_candidate, &old_candidate)) {
            //found a trivial match
            if (new_candidate.dist < old_candidate.dist) {
                //old candidate is trivial match (new is better), so disregard old and register the new one
                trivial_matches_map->push_candidate(new_candidate);
                std::cout << "Trivial match:" << std::endl;
                std::cout << "    old loc: " << std::setw(10) << old_candidate.loc << ", old query loc: " << std::setw(10) << old_candidate.query_loc << ", dist " << std::setw(10) << old_candidate.dist << std::endl;
                std::cout << "    new loc: " << std::setw(10) << new_candidate.loc << ", new query loc: " << std::setw(10) << new_candidate.query_loc << ", dist " << std::setw(10) << new_candidate.dist << std::endl;
            }
        } else {
            //no trivial match
            trivial_matches_map->push_candidate(new_candidate);
        }
    }

    void run() {
        TrivialMatchMap trivial_match_map;
        num_trivial_matches = 0;
        num_all_matches = 0;
        while ( ! input_pq.empty()) {
            InputFile* input_file = input_pq.top();
            input_pq.pop();
            size_t curr_query_pos = input_file->get_curr_query_pos();
            size_t prev_query_pos;
            do {
                //std::cout << "processing query " << curr_query_pos << ", candidate " << input_file->get_curr_candidate_pos() << ", dist " << input_file->get_curr_distance() << std::endl;
                process_input_file_line(*input_file, &trivial_match_map);
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
