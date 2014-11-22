#include "postprocessor.h"

boost::regex PostProcessor::InputFile::query_cand_dist_expr = boost::regex(R"(^(\d+)[ ,](\d+)[ ,](\d*\.?\d*))");
