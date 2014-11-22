#ifndef _UCR_H_
#define _UCR_H_

/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected © 2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/


#include "common.h"
#include "subsequence.h"

#define dist(x,y) ((x-y)*(x-y))

#define INF 1e20       //Pseudo Infitinte number for this code

class UCR_DTW {
private:
    SubsequenceLookup& subsequences;

    /// Data structure for sorting the query
    typedef struct Index
        {   double value;
            int    index;
        } Index;

    /// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
    static int comp(const void *a, const void* b)
    {   Index* x = (Index*)a;
        Index* y = (Index*)b;
        return fabs(y->value) - fabs(x->value);   // high to low
    }

    /// Calculate quick lower bound
    /// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
    /// However, because of z-normalization the top and bottom cannot give siginifant benefits.
    /// And using the first and last points can be computed in constant time.
    /// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
    double lb_kim_hierarchy(Subsequence const& query, Subsequence const& candidate, int len, double bsf = INF) const {
        /// 1 point at front and back
        double d, lb;
        double x0 = candidate.series_normalized[0];
        double y0 = candidate.series_normalized[len - 1];
        lb = dist(x0, query.series_normalized[0]) + dist(y0, query.series_normalized[len - 1]);
        if (lb >= bsf)   return lb;

        /// 2 points at front
        double x1 = candidate.series_normalized[1];
        d = std::min(dist(x1, query.series_normalized[0]), dist(x0, query.series_normalized[1]));
        d = std::min(d, dist(x1, query.series_normalized[1]));
        lb += d;
        if (lb >= bsf)   return lb;

        /// 2 points at back
        double y1 = candidate.series_normalized[len - 2];
        d = std::min(dist(y1, query.series_normalized[len - 1]), dist(y0, query.series_normalized[len - 2]));
        d = std::min(d, dist(y1, query.series_normalized[len - 2]));

        lb += d;
        if (lb >= bsf)   return lb;

        /// 3 points at front
        double x2 = candidate.series_normalized[2];
        d = std::min(dist(x0, query.series_normalized[2]), dist(x1, query.series_normalized[2]));
        d = std::min(d, dist(x2, query.series_normalized[2]));
        d = std::min(d, dist(x2, query.series_normalized[1]));
        d = std::min(d, dist(x2, query.series_normalized[0]));
        lb += d;
        if (lb >= bsf)   return lb;

        /// 3 points at back
        double y2 = candidate.series_normalized[len - 3];
        //XXX I think the following line is wrong: it doesn't match the previous patterns
        d = std::min(dist(y0, query.series_normalized[len - 3]), dist(y1, query.series_normalized[len - 3]));
        d = std::min(d, dist(y2, query.series_normalized[len - 3]));
        d = std::min(d, dist(y2, query.series_normalized[len - 2]));
        d = std::min(d, dist(y2, query.series_normalized[len - 1]));
        lb += d;

        return lb;
    }

    /// LB_Keogh 1: Create Envelop for the query
    /// Note that because the query is known, envelop can be created once at the begenining.
    ///
    /// Variable Explanation,
    /// order : sorted indices for the query.
    /// uo, lo: upper and lower envelops for the query, which already sorted.
    /// t     : a circular array keeping the current data.
    /// j     : index of the starting location in t
    /// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
    double lb_keogh_cumulative(Subsequence const& candidate, int* order, double *uo, double *lo, double *cb, int len, double best_so_far)
    {
        double lb = 0;
        double x, d;

        for (int i = 0; i < len && lb < best_so_far; i++)
        {
            x = candidate.series_normalized[order[i]];
            d = 0;
            if (x > uo[i])
                d = dist(x,uo[i]);
            else if(x < lo[i])
                d = dist(x,lo[i]);
            lb += d;
            cb[order[i]] = d;
        }
        return lb;
    }

    /// LB_Keogh 2: Create Envelop for the data
    /// Note that the envelops have been created (in main function) when each data point has been read.
    ///
    /// Variable Explanation,
    /// qo: sorted query
    /// cb: (output) current bound at each position. Used later for early abandoning in DTW.
    /// l,u: lower and upper envelop of the current data
    double lb_keogh_data_cumulative(int* order, double *qo, double *cb, const double *l, const double *u, int len, double mean, double stddev, double best_so_far = INF)
    {
        double lb = 0;
        double uu,ll,d;

        for (int i = 0; i < len && lb < best_so_far; i++)
        {
            uu = (u[order[i]]-mean)/stddev;
            ll = (l[order[i]]-mean)/stddev;
            d = 0;
            if (qo[i] > uu)
                d = dist(qo[i], uu);
            else
            {   if(qo[i] < ll)
                d = dist(qo[i], ll);
            }
            lb += d;
            cb[order[i]] = d;
        }
        return lb;
    }

    /// Print function for debugging
    void printArray(double *x, int len)
    {   for(int i=0; i<len; i++)
            printf(" %6.2lf",x[i]);
        printf("\n");
    }

    /// If expected error happens, teminated the program.
    void error(int id)
    {
        if(id==1)
            printf("ERROR : Memory can't be allocated!!!\n\n");
        else if ( id == 2 )
            printf("ERROR : File not Found!!!\n\n");
        else if ( id == 3 )
            printf("ERROR : Can't create Output File!!!\n\n");
        exit(1);
    }
public:
    UCR_DTW(SubsequenceLookup& subsequence_lookup): subsequences(subsequence_lookup) {}

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

    void single_pass(size_t query_position,
                     size_t candidate_increment,
                     std::function<dist_type()> weakest_distance_callback,
                     std::function<void(size_t, dist_type)> register_candidate_callback)
    {
        Subsequence const& query = subsequences[query_position];
        const unsigned int m = QUERY_LEN;
        if (query.range < MIN_RANGE) {
            msgl("Range is too small, skipping");
            return;
        }

        double bsf;          /// best-so-far
        //double *u, *l;
        double qo[m];
        double uo[m];
        double lo[m];
        int order[m];
        double cb[m];
        double cb1[m];
        double cb2[m];

        Index Q_tmp[m];



        long long i;
        double mean, stddev;
        int kim = 0,keogh = 0, keogh2 = 0;
        double dist=0, lb_kim=0, lb_k=0, lb_k2=0;

        bsf = INF;
        i = 0;
        const std::vector<double>& q = query.series_normalized;

        /// Create envelop of the query: lower envelop, l, and upper envelop, u
        const double* l = query.lemire_envelope.lower;
        const double* u = query.lemire_envelope.upper;

        /// Sort the query one time by abs(z-norm(q[i]))
        //todo: consider storing sorted query in cache (will this help anything?)
        for( i = 0; i<m; i++)
        {
            Q_tmp[i].value = q[i];
            Q_tmp[i].index = i;
        }
        qsort(Q_tmp, m, sizeof(Index), &comp);

        /// also create another arrays for keeping sorted envelop
        for( i=0; i<m; i++)
        {   int o = Q_tmp[i].index;
            order[i] = o;
            qo[i] = q[o];
            uo[i] = u[o];
            lo[i] = l[o];
        }

        /// Initial the cummulative lower bound
        for( i=0; i<m; i++)
        {   cb[i]=0;
            cb1[i]=0;
            cb2[i]=0;
        }

        int k=0;

        std::cerr << "starting: " << query_position << std::endl;
        for (size_t candidate_position = query_position + m; candidate_position < TIME_SERIES_LEN - m; candidate_position += candidate_increment) {
            Subsequence const& candidate = subsequences[candidate_position];
            if (candidate.range < MIN_RANGE) {
                continue;
            }
            const double* l_buff = candidate.lemire_envelope.lower;
            const double* u_buff = candidate.lemire_envelope.upper;
            mean = candidate.mean;
            stddev = candidate.stddev;

            /// Use a constant lower bound to prune the obvious subsequence
            lb_kim = lb_kim_hierarchy(query, candidate, m, bsf);


            if (lb_kim < bsf)
            {
                /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                /// uo, lo are envelop of the query.
                lb_k = lb_keogh_cumulative(candidate, order, uo, lo, cb1, m, bsf);
                if (lb_k < bsf)
                {
                    const std::vector<double>& tz = candidate.series_normalized;

                    /// Use another lb_keogh to prune
                    /// qo is the sorted query. tz is unsorted z_normalized data.
                    /// l_buff, u_buff are big envelop for all data in this chunk
                    lb_k2 = lb_keogh_data_cumulative(order, qo, cb2, l_buff, u_buff, m, mean, stddev, bsf);
                    if (lb_k2 < bsf)
                    {
                        /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                        /// Note that cb and cb2 will be cumulative summed here (possible optimization opportunity - move to dtw function?).
                        if (lb_k > lb_k2)
                        {
                            cb[m-1]=cb1[m-1];
                            for(k=m-2; k>=0; k--)
                                cb[k] = cb[k+1]+cb1[k];
                        }
                        else
                        {
                            cb[m-1]=cb2[m-1];
                            for(k=m-2; k>=0; k--)
                                cb[k] = cb[k+1]+cb2[k];
                        }


                        /// Compute DTW and early abandoning if possible
                        dist = dtw(tz, q, cb, m, bsf);

                        if( dist < bsf )
                        {   /// Update bsf
                            /// loc is the real starting location of the nearest neighbor in the file
                            register_candidate_callback(candidate_position, dist);
                            bsf = weakest_distance_callback();
                        }
                    } else
                        keogh2++;
                } else
                    keogh++;
            } else
                kim++;
        }
        std::cerr << "done" << std::endl;
        /*

        t2 = clock();
        printf("\n");

        /// Note that loc and i are long long.
        cout << "Location : " << loc << endl;
        cout << "Distance : " << sqrt(bsf) << endl;
        cout << "Data Scanned : " << i << endl;
        cout << "Total Execution Time : " << (t2-t1)/CLOCKS_PER_SEC << " sec" << endl;

        /// printf is just easier for formating ;)
        printf("\n");
        printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / i)*100);
        printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / i)*100);
        printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / i)*100);
        printf("DTW Calculation     : %6.2f%%\n", 100-(((double)kim+keogh+keogh2)/i*100));
        */
    }
};
#endif
