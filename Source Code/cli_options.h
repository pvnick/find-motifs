#ifndef _CLI_OPTIONS_H_
#define _CLI_OPTIONS_H_

#include <stdexcept>
#include <boost/program_options.hpp>

namespace program_options = boost::program_options;

class CLIOptions: public program_options::variables_map {
private:
    static CLIOptions instance;
    CLIOptions() = default;
public:
    static bool initialized;
    static void init(int argc, char *argv[]) {
        if ( ! initialized) {
            std::string cmd;
            // Declare the supported options.
            program_options::options_description desc("Allowed options");
            desc.add_options()
                ("K",
                    program_options::value<unsigned int>()->default_value(100)->value_name("INT"),
                    "set the number of top candidates to store per query")
                ("input-files",
                    program_options::value<std::vector<std::string>>()->multitoken()->value_name("[FILE...]")
                    , "list of results files for post-processing")
                ("command",
                    program_options::value<std::string>(&cmd)->required()->value_name("COMMAND"),
                    "command")
            ;

            try {
                program_options::store(program_options::command_line_parser(argc, argv)
                                            .options(desc)
                                            .run(),
                                       instance);

                program_options::notify(instance);
                if (cmd == "help")
                {
                    std::cout << desc << std::endl;
                    exit(0);
                }
            } catch(program_options::error& e) {
                std::cout << desc << std::endl;
                exit(0);
            }
            initialized = true;
        }
    }
    static CLIOptions const& get_instance() {
        if (initialized)
            return instance;
        throw std::logic_error("You must initialize the options class before retrieving an instance");
    }
};

#endif
