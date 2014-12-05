#ifndef _PIPELINE_H_
#define _PIPELINE_H_

#include <iostream>
#include <algorithm>
#include <functional>
#include <map>
#include <boost/any.hpp>

template<class Msg, class... Procs>
class Pipeline {
//recursion placeholder, do nothing
public:
    Pipeline(std::map<std::string, boost::any>& init_args) {}
    void go_not_first(Msg const& msg) {}
    void end() {}
};

/*
    A dead-simple data pipeline for chaining together state-maintaining processes in a producer-filter-consumer pattern.
    Supports one-to-zero, one-to-one, and one-to-many chaining.
    Definition is Pipeline<Msg, Proc, Procs...> pipeline(std::map<std::string, boost::any> init_args). Msg is the type of message
    being produced and passed to each Proc. Each Proc is either the producer (first proc in the parameter pack) a filter
    (second through second-to-last proc in the parameter pack), or the consumer (last proc in the parameter pack).
    init_args is a map of initialization arguments with string keys and boost::any values. init_args is passed as the only
    parameter to each Proc when it is constructed.

    Each Proc is a class that should expose, at a minimum, these public members:
        Proc(std::map<std::string, boost::any>& init_args)
            Called by the pipeline when constructing each proc, using the "new" keyword.
        void go(Msg const& msg, std::function<void(Msg const&)> yield)
            ...if the Proc is a filter/consumer (ie *not* the first Proc in the parameter pack).
            msg is what's received from the previous proc.
            yield passes a message to the next proc (typically unused by the consumer)
        void init(std::function<void(Msg const&)> yield)
            ...if the Proc is the producer (is the first Proc in the parameter pack)
            The producer is in charge of generating and passing each of the seed messages, typically from the set of
            initialization arguments.
        void end(std::function<void()> yield)
            Called when there are no more messages being generated, and is used for cleanup (eg filtering/consuming
            any left-over queued messages). Called after the producer's init() function returns. Each Proc
            should in turn call yield() to send the termination signal to the next process.

    Usage is pipeline.go(), which calls producer.init(yield). That function is expected to call yield()
    with each message it wishes to pass through the pipeline. After producer.init finishes, the pipeline calls
    producer.end(yield) to instruct it to pass the end signal up the chain.

*/
template<class Msg, class Proc, class... Procs>
class Pipeline<Msg, Proc, Procs...>: Pipeline<Msg, Procs...> {
private:
    Proc* processor;
protected:
    void go_not_first(Msg const& msg) {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Msg, Proc, Procs...>::sink, this, _1);
        processor->go(msg, sink);
    }
    void end() {
        using namespace std::placeholders;
        auto end_sink = std::bind(&Pipeline<Msg, Proc, Procs...>::end_sink, this);
        processor->end(end_sink);
    }
public:
    Pipeline(std::map<std::string, boost::any>& init_args): Pipeline<Msg, Procs...>(init_args), processor(new Proc(init_args)) {}
    void go() {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Msg, Proc, Procs...>::sink, this, _1);
        processor->init(sink);
        end();
    }
    ~Pipeline() {
        delete processor;
    }
    void sink(Msg const& msg) {
        Pipeline<Msg, Procs...>::go_not_first(msg);
    }
    void end_sink() {
        Pipeline<Msg, Procs...>::end();
    }
};

#endif
