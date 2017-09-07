#include "data_ioutil.h"

#include "../settings/settings.h"
#include "../data/data_description.h"

#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

const size_t DataIOUtil::AUDIO_FEATUREVECTOR_ENTRY_SIZE = 32;
const size_t DataIOUtil::AUDIO_FEATUREVECTOR_LENGTH = 42;
const std::string DataIOUtil::DATA_FILE_EXTENSION = ".csv";

void DataIOUtil::printStatistis(commonutil::DataSet const & dataset)
{
 	if(commonSettings().printStatistics)
	{
		idx_type d = dataset.d();
		idx_type n = dataset.n();
		
		fp_type min = (dataset.points.col(1)-dataset.points.col(2)).norm();
		fp_type max = min;
		Vector min_d = (dataset.points.col(1)-dataset.points.col(2)).cwiseAbs();
		Vector max_d = min_d;
		fp_type sum = 0;
		Matrix covar = Matrix::Zero(d,d);
		Vector mean = Vector::Zero(d);
		for(idx_type i=0; i<n; ++i)
			mean += dataset.points.col(i);
		mean /= n;
		for(idx_type i=0; i<n; ++i) // up to i=n due to covar computation!
		{
			Vector pi = dataset.points.col(i);
			Vector diff_pi_mean = pi-mean;
			covar += diff_pi_mean * diff_pi_mean.transpose();
			for(idx_type j=i+1; j<n; ++j)
			{
				Vector pj = dataset.points.col(j);
				for(idx_type l=0; l<d; ++l)
				{
					Vector diff = pi - pj;
					fp_type diffnorm = diff.norm();
					sum += diff.squaredNorm();
					min = fmin(min, diffnorm);
					max = fmax(max, diffnorm);
					min_d = min_d.cwiseMin(diff);
					max_d = max_d.cwiseMax(diff);
				}
			}
		}
		covar /= n;
		
		std::cout << "statistics of generated data set:" 
		          << "  min_{n,m} ||x_n-x_m|| = " << min << std::endl
			  << "  (min_{n,m} |x_n-x_m|_d)_d = " << min_d.transpose() << std::endl
			  << "  (max_{n,m} |x_n-x_m|_d)_d = " << max_d.transpose() << std::endl
			  << "  *Spread(X)* = max_{n,m} ||x_n-x_m|| = " << max << std::endl
			  << "  *sigma^2(X)* = sum/(2N^2) = " << sum / (2*n*n) << std::endl
			  << "  *fnorm(Sigma^2(X))* =  = " << covar.norm() << std::endl;
			  
			  
	}
}

void DataIOUtil::runDataGeneration(std::string const& dataDescriptionFile, std::string const& outputDir)
{
	std::vector<DataDescription*> dataDescriptors = DataDescription::createDataDescriptions(dataDescriptionFile);

	IOUtil outil(".csv");
	std::stringstream pathStream;
	if(!outputDir.empty())
	{
		pathStream << outputDir << "/";
		boost::filesystem::path outPath(pathStream.str());
		boost::filesystem::create_directories(outPath);
	}
	
	commonutil::DataSet input;
	

	// loop over all data sets
	for (unsigned int dataIndex=0; dataIndex<dataDescriptors.size(); ++dataIndex)
	{
		// clear input
		input.weights.resize(0);
		input.points.resize(0,0);
	
		std::cout << ">" << std::endl
				  << ">> retrieving new data set from descriptor: "
				  << dataDescriptors[dataIndex]->tag() << std::endl
				  << ">" << std::endl;

		// load or generate data from descriptor		
		dataDescriptors[dataIndex]->retrieve(input);
		if (input.weights.size()==0)
		{
			std::cout << "   skipping empty data set." << std::endl << std::endl;
			continue;
		}
		
		idx_type n = input.points.cols();
		idx_type d = input.points.rows();
		
		
		if(!outputDir.empty())
		{
			outil.open(pathStream.str() + dataDescriptors[dataIndex]->tag());
			for(idx_type i=0; i<n; i++)
				outil.store(input.points.col(i));
			outil.close();
			
			std::cout << ">> stored " << pathStream.str() + dataDescriptors[dataIndex]->tag() << " in " << outputDir << std::endl << ">" << std::endl;
		}
		
		printStatistis(input);
	}
}



void DataIOUtil::loadDataFromCSV(commonutil::DataSet& target, std::string const& filename,
	unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize)
{
	std::ifstream filestream;

	// in the first pass we count the number of lines and check, if all input points have the
	// same dimension or if mindim > 0, if all have at least dimension mindim
	
	std::string line;
	idx_type num = 0;
	idx_type dim = numCoords==0?0:firstCoord+numCoords-1;
	if (commonSettings().verbose)
		std::cout << "analyzing file " << filename << std::endl;
	filestream.open(filename.c_str());
	while(getline(filestream, line))
	{
		num++;
		std::istringstream iss(line);
		idx_type tmpdim = 0;
		std::string entry;
		while (getline(iss, entry, IOUtil::SEPARATOR))
		{
			++tmpdim;
			if (tmpdim==dim)
				break;
		}


		if (dim != tmpdim && dim != 0)
		{
			std::cout << "ioutil::loadDataFromCSV - found data points with different dimensions "
					  << dim << " and " << tmpdim << std::endl;
			std::cout << "in line number " << num << "... " << line << std::endl;
			filestream.close();
			gmmlab_throw("ioutil::loadDataFromCSV - loading failed!!!");
		}
	}
	filestream.close();

	if (firstCoord>dim)
	{
		if (commonSettings().verbose)
			std::cout << "ioutil::loadDataFromCSV - selected coordinates out of range!" << std::endl;
		return;
	}

	// dimension
	unsigned int coordsToBeLoaded = numCoords;
	if (numCoords==0 || firstCoord+numCoords-1>dim)
		coordsToBeLoaded = dim-firstCoord+1;

	idx_type shiftlen = (blockLength-1)*coordsToBeLoaded;

	unsigned int targetDim = blockLength*coordsToBeLoaded;
	if (target.weights.size()==0)
		target.points.resize(targetDim, 0);
	if (target.points.rows()!=targetDim)
		gmmlab_throw("target dataset has incompatible dimension.");

	// resize the input matrix so that each data point will form one column
	unsigned int blocks = num-blockLength+1;
	idx_type oldsize = target.weights.size();
	idx_type newsize = target.weights.size()+blocks;
	target.weights.conservativeResize(newsize);
	target.weights.tail(blocks) = Vector::Ones(blocks);
	target.points.conservativeResize(targetDim, newsize);

	// in the second pass we read the coordinates blockwise into the input matrix
	
	std::cout << "reading " << blocks << " data points in " << targetDim << " dimensions from "
			  << num << " lines of " << filename << std::endl;
	Vector oneblock(targetDim);
	filestream.open(filename.c_str());
	for (idx_type i=0; i<num; ++i)
	{
		getline(filestream, line);
		std::istringstream iss(line);
		std::string entry;
		oneblock.head(shiftlen) = oneblock.tail(shiftlen);
		for (idx_type j=0; j<dim; ++j)
		{
			getline(iss, entry, ',');
			if(j>=firstCoord-1 && j<firstCoord-1+coordsToBeLoaded)
				oneblock(targetDim-coordsToBeLoaded+j-firstCoord+1)=atof(entry.c_str());
		}
		if (i>=blockLength-1)
			target.points.col(oldsize+i) = oneblock;
//		std::cout << "new block: " << oneblock << std::endl;
	}
	filestream.close();

	// normalization
	if (normalize)
	{
		std::cout << "   normalizing data..." << std::endl;
		for (idx_type i=0; i<targetDim; ++i)
		{
			fp_type coeff = target.points.row(i).tail(blocks).minCoeff();
			target.points.row(i).tail(blocks) -= Vector::Constant(blocks, coeff);
			coeff = target.points.row(i).tail(blocks).maxCoeff();
			if (coeff>0)
				target.points.row(i).tail(blocks) /= coeff;
		}
	}
//	std::cout << std::endl;
}

// loads fe file and returns number of dropped frames at the beginning
idx_type DataIOUtil::loadDataFromFE(commonutil::DataSet& target, std::string const& filename,
	unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize, bool noFrameDrop)
{
	if (firstCoord>=AUDIO_FEATUREVECTOR_LENGTH)
	{
		if (commonSettings().verbose)
			std::cout << "ioutil::loadDataFromFE - selected coordinates out of range!" << std::endl;
		return 0;
	}

	// dimension
	unsigned int coordsToBeLoaded = numCoords;
	if (numCoords==0 || firstCoord+numCoords>AUDIO_FEATUREVECTOR_LENGTH)
		coordsToBeLoaded = AUDIO_FEATUREVECTOR_LENGTH-firstCoord;
	std::streamoff pregap = firstCoord*sizeof(float);
	std::streamoff postgap = (AUDIO_FEATUREVECTOR_LENGTH-firstCoord-coordsToBeLoaded)*sizeof(float);

	idx_type shiftlen = (blockLength-1)*coordsToBeLoaded;
	
	unsigned int targetDim = blockLength*coordsToBeLoaded;
	if (target.weights.size()==0)
		target.points.resize(targetDim, 0);
	if (target.points.rows()!=targetDim)
		gmmlab_throw("target dataset has incompatible dimension.");

	// open file
	if(commonSettings().verbose)
		std::cout << "   reading from file " << filename << std::endl;
	std::ifstream ifs(filename.c_str(), std::ifstream::binary | std::ifstream::in);
	if(!ifs)
		gmmlab_throw("Failed to open file.");

	// determine number of entries
	ifs.seekg(0, std::ifstream::end);
	idx_type numberOfEntries = ifs.tellg()/sizeof(float);
	if (numberOfEntries % AUDIO_FEATUREVECTOR_LENGTH)
		gmmlab_throw("Wrong format.");
	idx_type numberOfVectors = numberOfEntries/AUDIO_FEATUREVECTOR_LENGTH;
	ifs.seekg(0, std::ifstream::beg);
	if(commonSettings().verbose)
		std::cout << "      " << numberOfVectors << " feature vectors found." << std::endl;

	Matrix matrix(targetDim, numberOfVectors);
	float buffer;

	// read entries
	Vector first;
	idx_type dropped = 0;
	bool different = false; // indicates whether there has been a point that differs completely from the first point
	if(noFrameDrop)
	  different = true;
	Vector block(targetDim);
	unsigned int blocks = 0;
	for(std::size_t i=0; i < numberOfVectors; ++i)
	{
		// read the next point
		Vector next(coordsToBeLoaded);
		
		ifs.seekg(pregap, std::ios_base::cur);
		for(std::size_t j=0; j<coordsToBeLoaded; ++j)
		{
			ifs.read((char*) &buffer, sizeof(float));
			next(j) = buffer;
		}
		ifs.seekg(postgap, std::ios_base::cur);

		// datawash
		if(i==0)
			first = next;
		else
		{
			if(!different && commonutil::componentWiseDifferent(first, next))
			{
				different = true;
				dropped = i;
				if(commonSettings().verbose)
					std::cout << "      " << dropped << " vectors skipped." << std::endl;
			}

			if(different)
			{
				block.head(shiftlen) = block.tail(shiftlen);
				block.tail(coordsToBeLoaded) = next;
				if (i-dropped >= blockLength-1)
					matrix.col(i-dropped-blockLength+1) = block;
			}
		}
	}

	// close file
	ifs.close();


	idx_type numread = numberOfVectors-dropped-blockLength+1;
	matrix.conservativeResize(targetDim, numread);

	// normalization
	if (normalize)
	{
		std::cout << "   normalizing data..." << std::endl;
		for (idx_type i=0; i<targetDim; ++i)
		{
			fp_type coeff = matrix.row(i).minCoeff();
			matrix.row(i) -= Vector::Constant(numread, coeff);
			coeff = matrix.row(i).maxCoeff();
			if (coeff>0)
				matrix.row(i) /= coeff;
		}
	}
//	std::cout << std::endl;

	idx_type newsize = target.weights.size()+numread;
	target.weights.conservativeResize(newsize);
	target.weights.tail(numread) = Vector::Ones(numread);
	target.points.conservativeResize(targetDim, newsize);
	target.points.rightCols(numread) = matrix;

	return dropped;
}

void DataIOUtil::loadData(commonutil::DataSet& target, std::string const& pattern,
	unsigned int firstCoord, unsigned int numCoords, unsigned int blockLength, bool normalize, bool noFrameDrop)
{
	std::cout << "Loading datapoints from files starting with: " << pattern << std::endl;
	idx_type presize = target.weights.size();

	// get filenames from pattern and sort alphabetically
	std::vector<std::string> filelist;
	IOUtil::getFilenames(pattern, filelist);
	std::sort(filelist.begin(), filelist.end());

	if (filelist.empty())
		std::cout << "   WARNING - No files found!" << std::endl;
	else
	{
		// check for common file extension to apply proper load procedure
		std::string common = IOUtil::commonFileExtension(filelist);
		boost::to_upper(common);
		IOUtil::FileFormat ff = IOUtil::FileFormat::UNKNOWN;
		if (common == "CSV")
			ff = IOUtil::FileFormat::CSV;
		else if (common == "LAB")
			ff = IOUtil::FileFormat::LAB;
		else if (common == "FE")
			ff = IOUtil::FileFormat::FE;
		else
			std::cout << "   WARNING - Unknown or diverse file extension " << common << " in filelist of size "<< filelist.size() << "!" << std::endl;
		
		// iterate over list and load files
		idx_type maxdrop = 0;
		for (std::vector<std::string>::iterator it = filelist.begin(); it!=filelist.end(); ++it)
			switch (ff)
			{
				case IOUtil::FileFormat::CSV:
					loadDataFromCSV(target, *it, firstCoord, numCoords, blockLength, normalize);
					break;
// 				case IOUtil::FileFormat::LAB:
// 					loadDataFromLAB(target, *it, firstCoord, numCoords, blockLength, normalize);
// 					break;
				case IOUtil::FileFormat::FE:
					idx_type drop = loadDataFromFE(target, *it, firstCoord, numCoords, blockLength, normalize, noFrameDrop);
					if (drop>maxdrop)
						maxdrop = drop;
					break;
			}
		if (ff==IOUtil::FileFormat::FE)
			std::cout << "   maximum frame drop is " << maxdrop << "." << std::endl;
		std::cout << "   " << target.weights.size()-presize << " data points read." << std::endl;
	}
}

void DataIOUtil::store ( const Vector& entry )
{
	IOUtil::store(entry);
}

void DataIOUtil::store ( commonutil::DataSet const& input )
{
	for(idx_type i = 0; i < input.n(); ++i)
		store(input.points.col(i));
}