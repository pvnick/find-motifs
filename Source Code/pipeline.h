#ifndef _PIPELINE_H_
#define _PIPELINE_H_

#include <iostream>
#include <algorithm>
#include <functional>
#include <map>
#include <boost/any.hpp>

//todo: push this to its own repo in github in case others find it useful.
//although, i feel like it's a reinvention of some wheel, somewhere, that i was unable to find on google.

template<class... Procs>
class Pipeline {
//recursion placeholder, do nothing
public:
    Pipeline(std::map<std::string, boost::any> init_args) {}
    void go_not_first(boost::any const& msg) {}
    void end() {}
};

/*
    A dead-simple data pipeline for chaining together state-maintaining processes in a producer-filter-consumer pattern.
    Supports one-to-zero, one-to-one, and one-to-many chaining.
    Definition is Pipeline<Proc, Procs...> pipeline(std::map<std::string, boost::any> init_args). Each Proc is either the
    producer (first proc in the parameter pack) a filter (second through second-to-last proc in the parameter pack), or
    the consumer (last proc in the parameter pack).
    init_args is a map of initialization arguments with string keys and boost::any values. init_args is passed as the only
    parameter to each Proc when it is constructed.

    Each Proc is a class that should expose, at a minimum, these public members:
        Proc(std::map<std::string, boost::any> init_args)
            Called by the pipeline when constructing each proc, using the "new" keyword.
        void go(boost::any const& msg, std::function<void(boost::any const&)> yield)
            ...if the Proc is a filter/consumer (ie *not* the first Proc in the parameter pack).
            msg is what's received from the previous proc.
            yield passes a message to the next proc (typically unused by the consumer)
        void init(std::function<void(boost::any const&)> yield)
            ...if the Proc is the producer (is the first Proc in the parameter pack).
            The producer is in charge of generating and passing each of the original messages, typically dependent upon the
            set of initialization arguments.
        void end(std::function<void(boost::any const&)> yield, std::function<void()> forward_end_signal)
            Called after the producer's init() function returns, indicating no more messages are coming from the producer.
            Used for cleanup (eg filtering/consuming any left-over queued messages). Each proc may yield messages
            to pass through the pipeline as part of their cleanup procedure. Afterwards, they should call forward_end_signal()
            to pass the termination signal to the next process.

    Usage is pipeline.go(), which calls producer.init(yield). That function is expected to call yield()
    with each message it wishes to pass through the pipeline. After producer.init finishes, the pipeline calls
    producer.end(yield, forward_end_signal) to instruct it to pass the end signal up the chain.

*/
template<class Proc, class... Procs>
class Pipeline<Proc, Procs...>: Pipeline<Procs...> {
private:
    Proc* processor;
protected:
    void go_not_first(boost::any const& msg) {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Proc, Procs...>::sink, this, _1);
        processor->go(msg, sink);
    }
    void end() {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Proc, Procs...>::sink, this, _1);
        auto forward_end_signal = std::bind(&Pipeline<Proc, Procs...>::end_sink, this);
        processor->end(sink, forward_end_signal);
    }
public:
    Pipeline(std::map<std::string, boost::any> init_args): Pipeline<Procs...>(init_args), processor(new Proc(init_args)) {}
    void go() {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Proc, Procs...>::sink, this, _1);
        processor->init(sink);
        end();
    }
    ~Pipeline() {
        delete processor;
    }
    void sink(boost::any const& msg) {
        Pipeline<Procs...>::go_not_first(msg);
    }
    void end_sink() {
        Pipeline<Procs...>::end();
    }
};

#endif
