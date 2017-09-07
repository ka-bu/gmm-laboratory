#ifndef IOUTIL_H
#define IOUTIL_H

#include "../base.h"
#include <fstream>

#include <vector>

/**
 * @brief 
 */
class IOUtil
{
public:
	
	enum FileFormat { UNKNOWN, CSV, LAB, FE, CLS, GMM };
	
	static const char SEPARATOR;
	static const char COMMENT;
	static const std::string LIST_FILE_EXTENSION;
	static const std::string NA;
		
	static void getFilenames(std::string const& pattern, std::vector<std::string>& list);
	static const std::string commonFileExtension(std::vector<std::string> const& list);
	static std::vector<std::string> readDirectory(const std::string&);
	static void set_precision(std::fstream& filestream);
	static void set_precision(std::ofstream& filestream);
	static bool ends_with(std::string s, std::string postfix);
	
	
	IOUtil();
	IOUtil(std::string file_extension);
	~IOUtil();
	void open(std::string const& path, unsigned int k = 0);
	bool is_open();
	virtual void header();
	virtual void store(Vector const& entry);
	void close();
	
protected:
	
	std::fstream filestream;
	
	unsigned int k;
	std::string file_extension = ".csv";
	
};

#endif