#include "gen_umrc_description.h"

#include "../settings/configparser.h"
#include "../settings/settings.h"
#include "datageneration.h"

GenUMRCDescription::GenUMRCDescription(unsigned int n, unsigned int d, unsigned int k, fp_type weightFloatExp, fp_type separation, fp_type sqrtEWFloatExp, fp_type minMinSqrtEW, fp_type maxMinSqrtEW, fp_type minSqrtEWProp, fp_type maxSqrtEWProp, std::string measurementError, uint32_t seed) :
size(n), dim(d), comp(k), weightFloatExp(weightFloatExp), separation(separation), sqrtEWFloatExp(sqrtEWFloatExp), minMinSqrtEW(minMinSqrtEW), maxMinSqrtEW(maxMinSqrtEW), minSqrtEWProp(minSqrtEWProp), maxSqrtEWProp(maxSqrtEWProp), measurementError(measurementError), genSeed(seed)
{
	std::stringstream sstream;
	sstream << "gen_umrc_n" << n << "_d" << d << "_k" << k  << "_sep" << separation << "_w" << weightFloatExp << "_ewexp" << sqrtEWFloatExp << "_ewmin" << minMinSqrtEW << "_ewmax" << maxMinSqrtEW << "_ewpmin" << minSqrtEWProp << "_ewpmax" << maxSqrtEWProp << "_me" << measurementError << "_g" << seed;
	this->nametag = sstream.str();
	
	this->kList.push_back(k);
}

Parameters GenUMRCDescription::retrieve(commonutil::DataSet& input)
{
	std::mt19937 gen(this->genSeed);
	
	if(commonSettings().verbose)
		std::cout << "generate random gmm according to... " << std::endl
		          << "   description = " << this->nametag << std::endl;
	
// 	Parameters truth = datageneration::generateRandomMPEMMWithUniformMeans(dim, comp, gen, 1, 1,
// 									       separation,
// 													   minMinSqrtEW, maxMinSqrtEW, 
// 													   minSqrtEWProp, maxSqrtEWProp,
// 													   weightFloatExp, sqrtEWFloatExp);
	Parameters truth = datageneration::generateRandomGMMWithUniformMeans(dim, comp, gen, 
									       separation,
													   minMinSqrtEW, maxMinSqrtEW, 
													   minSqrtEWProp, maxSqrtEWProp,
													   weightFloatExp, sqrtEWFloatExp);
	unsigned int uniformNoisePoints = round(0.1 * this->size);
	
	if(this->measurementError ==  datageneration::ME_NO_ME) 
		datageneration::generateInputFromGMM(input, truth, this->size, gen);
	else if(this->measurementError == datageneration::ME_BLUR_CORRELATED)
			datageneration::generateCorrelatedNoisyInputFromGMM(input, truth, this->size, 0.5, gen);
	else if(this->measurementError == datageneration::ME_BLUR_DIFFICULT)
			datageneration::generateDifficultNoisyInputFromGMM(input, truth, this->size, gen);
	else if(this->measurementError == datageneration::ME_UNIFORM_NOISE)
	{
			datageneration::generateInputFromGMM(input, truth, this->size - uniformNoisePoints, gen);
			datageneration::addUniformNoise(input, uniformNoisePoints, gen);
	}
	else if(this->measurementError == datageneration::ME_GAUSSIAN_NOISE)
	{
		gmmlab_throw(" To be implemented");
	}
	else
	{
		std::cout << "ERROR:  Measurement error " << this->measurementError << std::endl;
		gmmlab_throw(" Unknown type of measurement error");
	}

	return truth;
}



