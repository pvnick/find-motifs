#ifndef _POSTPROCESS_H_
#define _POSTPROCESS_H_

#include <fstream>
#include <queue>
#include <boost/regex.hpp>

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
        bool is_alive;
        static boost::regex query_cand_dist_expr; //matches input file line fields (UINT UINT DOUBLE)
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
            while (is_alive) {
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
                    is_alive = false; //end of input file
            }
        }
        InputFile(std::string const& input_file_name): input(input_file_name), is_alive(true) {
            next();
        }
        size_t get_curr_query_pos() { return curr_query_pos; }
        size_t get_curr_candidate_pos() { return curr_candidate_pos; }
        double get_curr_distance() { return curr_distance; }
    };

    PointerPriorityQueue<InputFile> input_pq;
public:
    PostProcessor(std::vector<std::string> const& input_file_names) {
        for (std::string const& input_file_name: input_file_names) {
            input_pq.push(new InputFile(input_file_name));
        }
    }
    void run() {
        while ( ! input_pq.empty()) {
            InputFile* input_file = input_pq.top();
            size_t curr_query_pos = input_file->get_curr_query_pos();
            size_t prev_query_pos;
            do {
                std::cout << "processing query " << curr_query_pos << ", candidate " << input_file->get_curr_candidate_pos() << ", dist " << input_file->get_curr_distance() << std::endl;
                input_file->next();
                prev_query_pos = curr_query_pos;
                curr_query_pos = input_file->get_curr_query_pos();
            } while (prev_query_pos == curr_query_pos);
            input_pq.pop();
            input_pq.push(input_file);
        }
    }
};

#endif
