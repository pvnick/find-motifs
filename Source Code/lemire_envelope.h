#ifndef _LEMIRE_ENVELOPE_H_
#define _LEMIRE_ENVELOPE_H_

#include "common.h"
class LemireEnvelope {
private:
    /// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
    struct deque
    {   int *dq;
        int size,capacity;
        int f,r;
    };

    /// Initial the queue at the begining step of envelop calculation
    void init(deque *d, int capacity)
    {
        d->capacity = capacity;
        d->size = 0;
        d->dq = new int[d->capacity];
        d->f = 0;
        d->r = d->capacity-1;
    }

    /// Destroy the queue
    void destroy(deque *d)
    {
        delete[] d->dq;
    }

    /// Insert to the queue at the back
    void push_back(struct deque *d, int v)
    {
        d->dq[d->r] = v;
        d->r--;
        if (d->r < 0)
            d->r = d->capacity-1;
        d->size++;
    }

    /// Delete the current (front) element from queue
    void pop_front(struct deque *d)
    {
        d->f--;
        if (d->f < 0)
            d->f = d->capacity-1;
        d->size--;
    }

    /// Delete the last element from queue
    void pop_back(struct deque *d)
    {
        d->r = (d->r+1)%d->capacity;
        d->size--;
    }

    /// Get the value at the current position of the circular queue
    int front(struct deque *d)
    {
        int aux = d->f - 1;

        if (aux < 0)
            aux = d->capacity-1;
        return d->dq[aux];
    }

    /// Get the value at the last position of the circular queueint back(struct deque *d)
    int back(struct deque *d)
    {
        int aux = (d->r+1)%d->capacity;
        return d->dq[aux];
    }

    /// Check whether or not the queue is empty
    int empty(struct deque *d)
    {
        return d->size == 0;
    }

    /// Finding the envelop of min and max value for LB_Keogh
    /// Implementation idea is intoruduced by Danial Lemire in his paper
    /// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
    void lower_upper_lemire(const double *t, int len, int r, double *l, double *u)
    {
        struct deque du, dl;

        init(&du, 2*r+2);
        init(&dl, 2*r+2);

        push_back(&du, 0);
        push_back(&dl, 0);

        for (int i = 1; i < len; i++)
        {
            if (i > r)
            {
                u[i-r-1] = t[front(&du)];
                l[i-r-1] = t[front(&dl)];
            }
            if (t[i] > t[i-1])
            {
                pop_back(&du);
                while (!empty(&du) && t[i] > t[back(&du)])
                    pop_back(&du);
            }
            else
            {
                pop_back(&dl);
                while (!empty(&dl) && t[i] < t[back(&dl)])
                    pop_back(&dl);
            }
            push_back(&du, i);
            push_back(&dl, i);
            if (i == 2 * r + 1 + front(&du))
                pop_front(&du);
            else if (i == 2 * r + 1 + front(&dl))
                pop_front(&dl);
        }
        for (int i = len; i < len+r+1; i++)
        {
            u[i-r-1] = t[front(&du)];
            l[i-r-1] = t[front(&dl)];
            if (i-front(&du) >= 2 * r + 1)
                pop_front(&du);
            if (i-front(&dl) >= 2 * r + 1)
                pop_front(&dl);
        }
        destroy(&du);
        destroy(&dl);
    }
public:
    double lower[QUERY_LEN];
    double upper[QUERY_LEN];
    LemireEnvelope() = default;
    LemireEnvelope(const std::vector<double>& time_series, size_t start_position, int r) {
        lower_upper_lemire(time_series.data() + start_position, QUERY_LEN, r, lower, upper);
    }
    ~LemireEnvelope() {}
};

#endif
