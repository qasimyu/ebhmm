// ***************************************************************************
// InputParser.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _INPUTPARSER_H
#define _INPUTPARSER_H

class InputParser {
	private:
		void parseArgs_prepareInputs(int argc, char *argv[]);
		void parseArgs_callCNAs(int argc, char *argv[]);
		void usage(const char* app);
		void usage_prepareInputs(const char* app);
		void usage_callCNAs(const char* app);
	public:
		void parseArgs(int argc, char *argv[]);
};

#endif
