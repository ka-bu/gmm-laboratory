#ifndef MODEL_SELECTION_CRITERIA_H
#define MODEL_SELECTION_CRITERIA_H


#include "../../../base/commonutil.h"
#include "../../../base/parameters.h"

namespace model_selection
{

	fp_type bic(fp_type nll, idx_type n, idx_type k);
	fp_type bic(const commonutil::DataSet& dataset, const Parameters& gmmdesc, const std::vector<idx_type> & indices_of_points_to_be_used = std::vector<idx_type>());
	
	fp_type aic(fp_type nll, idx_type k);
	fp_type aic(const commonutil::DataSet& dataset, const Parameters& gmmdesc, const std::vector<idx_type> & indices_of_points_to_be_used = std::vector<idx_type>());
	
	fp_type mdl(fp_type nll, idx_type n, idx_type d, idx_type k);
	fp_type mdl(const commonutil::DataSet& dataset, const Parameters& gmmdesc, const std::vector<idx_type> & indices_of_points_to_be_used = std::vector<idx_type>());
	
	//fp_type minimum_description_length(commonutil::DataSet& dataset, GMMDesc& gmmdesc);
	
	
}



#endif