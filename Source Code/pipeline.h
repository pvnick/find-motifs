#ifndef _PIPELINE_H_
#define _PIPELINE_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

template<class Msg, class InitArgs, class... Procs>
class Pipeline {
//recursion placeholder, do nothing
public:
    Pipeline() {}
    void go_not_first(Msg const& msg) {}
    void end() {}
};

/*
    A dead-simple data pipeline for chaining together state-maintaining processes in a producer-filter-consumer pattern.
    Supports one-to-zero, one-to-one, and one-to-many chaining.
    Definition is Pipeline<Msg, InitArgs, Proc, Procs...> pipeline, where InitArgs is the type of initialization arguments
    passed to the producer to tell it how to produce the messages. Msg is the type of message being produced and
    passed to each Proc. The last Proc in the parameter pack is the consumer.

    Each Proc should expose these functions:
        void go(Msg const& msg, std::function<void(Msg const&)> yield)
            ...if the Proc is a filter/consumer (ie *not* the first Proc in the parameter pack).
            msg is what's received from the previous proc.
            yield passes a message to the next proc.
        void init(InitArgs const& args, std::function<void(Msg const&)> yield)
            ...if the Proc is the producer (is the first Proc in the parameter pack)
            args is the object (eg initializer_list/vector/etc) containing the list of arguments needed to produce messages
            yield passes a message to the next proc.
        void end(std::function<void()> yield)
            Called when there are no more messages being generated, and is used for cleanup (eg filtering/consuming
            any left-over queued messages). Called after the producer's init() function returns. Each Proc
            should in turn call yield() to send the termination signal to the next process.

    Since Proc instances are constructed using the "new" keyword without arguments, they should be
    constructable without arguments.

    Usage is pipeline.go(InitArgs args), which calls producer.init(args, yield). That function is expected to call yield()
    with each message it wishes to pass through the pipeline.
*/
template<class Msg, class InitArgs, class Proc, class... Procs>
class Pipeline<Msg, InitArgs, Proc, Procs...>: Pipeline<Msg, InitArgs, Procs...> {
private:
    Proc* processor;
protected:
    void go_not_first(Msg const& msg) {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Msg, InitArgs, Proc, Procs...>::sink, this, _1);
        processor->go(msg, sink);
    }
    void end() {
        using namespace std::placeholders;
        auto end_sink = std::bind(&Pipeline<Msg, InitArgs, Proc, Procs...>::end_sink, this);
        processor->end(end_sink);
    }
public:
    Pipeline(): Pipeline<Msg, InitArgs, Procs...>(), processor(new Proc) {}
    void go(InitArgs const& args) {
        using namespace std::placeholders;
        auto sink = std::bind(&Pipeline<Msg, InitArgs, Proc, Procs...>::sink, this, _1);
        processor->init(args, sink);
        end();
    }
    ~Pipeline() {
        delete processor;
    }
    void sink(Msg const& msg) {
        Pipeline<Msg, InitArgs, Procs...>::go_not_first(msg);
    }
    void end_sink() {
        Pipeline<Msg, InitArgs, Procs...>::end();
    }
};

#endif
