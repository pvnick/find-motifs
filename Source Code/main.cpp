//config
#define USE_MPI

//hardcoded values to prevent memory allocation
#define TIME_SERIES_LEN 1663231
#define QUERY_LEN 100
#define MIN_RANGE 5
#define WARPING_WINDOW 0.05
#define WARPING_r (WARPING_WINDOW <= 1) ? (const int)(WARPING_WINDOW * QUERY_LEN) : (const int)WARPING_WINDOW
#define SERIES_FILEPATH "/scratch/lfs/pvnick/oximetry.txt"

#include "common.h"
#include "find_motifs.h"
#include "ucr_dtw.h"

#ifdef USE_MPI
    bool use_mpi = true;
#else
    bool use_mpi = false;
#endif

size_t distributed_query_start_position(size_t series_length, unsigned int proc_rank, unsigned int num_procs) {
    size_t search_space_length = floor((float)(num_procs - proc_rank) / num_procs * series_length * (series_length + 1.0) / 2.0);
    size_t position = ceil(-1.0 / 2.0 * sqrt(8.0 * search_space_length + 1.0) + series_length + 1.0 / 2.0);
    return position;
}

unsigned int query_start_pos() {
    if (use_mpi) {
        mpi::communicator world;
        return distributed_query_start_position(TIME_SERIES_LEN, world.rank(), world.size());
    } else {
        return 0;
    }
}

unsigned int query_end_pos() {
    if (use_mpi) {
        mpi::communicator world;
        return distributed_query_start_position(TIME_SERIES_LEN, world.rank() + 1, world.size());
    } else {
        return TIME_SERIES_LEN - QUERY_LEN;
    }
}

std::string get_output_filename() {
    if (use_mpi) {
        mpi::communicator world;
        std::ostringstream filename;
        filename << "results" << world.rank() << ".out";
        return filename.str();
    } else {
        return "results.out";
    }
}

std::ostream& msg(std::string str) {
    if (use_mpi) {
        mpi::communicator world;
        std::cerr << "Proc " << world.rank() << ": " << str;
    } else {
        std::cerr << str;
    }
    return std::cerr;
}

std::ostream& msgl(std::string str) {
    return msg(str) << std::endl;
}


/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(const std::vector<double>& A, const std::vector<double>& B, double *cb, int m, double bsf = INF)
{
    const int r = WARPING_r;
    double buffer1[2 * r + 1];
    double buffer2[2 * r + 1];
    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    for(int k=0; k<2*r+1; k++)
        buffer1[k] = buffer2[k] = INF;

    double *cost = buffer1;
    double *cost_prev = buffer2;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;


    for (i=0; i<m; i++)
    {
        k = std::max(0,r-i);
        min_cost = INF;

        for(j=std::max(0,i-r); j<=std::min(m-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = std::min( std::min( x, y) , z) + dist(A[i],B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        //todo: consider removing DTW early abandoning
        if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
            return min_cost + cb[i+r+1];

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    return final_dtw;
}

#ifdef USE_PROFILER
    #include "profiler.h"
#endif


/// Main Function
int main(int argc, char *argv[])
{
#ifdef USE_MPI
    mpi::environment env(argc, argv);
#endif
    unsigned int K = 100;

    size_t start_pos = query_start_pos();
    size_t end_pos = query_end_pos();
#ifdef USE_PROFILER
    ProfilerStart("/tmp/profile");
#endif

    std::string results_file_path = std::string("/scratch/lfs/pvnick/motif_results/") + get_output_filename();
    std::string series_filepath = "/home/pvnick/oximetry/data/oximetry.txt";
    MotifFinder engine(series_filepath, TIME_SERIES_LEN, K, results_file_path, '\t');
    //engine.run(start_pos, end_pos, 1);
    std::vector<double> const& time_series = engine.get_timeseries();
    std::cout << time_series[2] << std::endl;
    SubsequenceLookup subsequences = engine.get_subsequence_lookup();
    Subsequence const& subsequence = subsequences[5];
    std::cout << subsequence.time_series_pos << std::endl;

#ifdef USE_PROFILER
    ProfilerStop();
#endif
    return 0;
}
