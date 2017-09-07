#include "model_selection.h"

#include "gmmutil.h"


fp_type model_selection::bic(fp_type nll, idx_type n, idx_type k)
{
	return 2*nll+k*(log(1.*n)); //+log(2*M_PI));
}


fp_type model_selection::aic(fp_type nll, idx_type k)
{
	return 2.*nll + 2.*k;
}

fp_type model_selection::mdl(fp_type nll, idx_type n, idx_type d, idx_type k)
{

	return 2.*nll+k*(d+d*(d+1.)/2.+1)*log(1.*n);
}

fp_type model_selection::bic(const commonutil::DataSet& dataset, const Parameters& gmmdesc, const std::vector< idx_type > & indices_of_points_to_be_used)
{
	unsigned int n = indices_of_points_to_be_used.empty() ? dataset.points.cols() : indices_of_points_to_be_used.size();
	unsigned int k = gmmdesc.weights.size();
	fp_type nll = gmmutil::nll(dataset, gmmdesc, indices_of_points_to_be_used);
	return model_selection::bic(nll, n, k);
}


fp_type model_selection::aic(const commonutil::DataSet& dataset, const Parameters& gmmdesc, const std::vector< idx_type > & indices_of_points_to_be_used)
{
	unsigned int k = gmmdesc.weights.size();
	fp_type nll = gmmutil::nll(dataset, gmmdesc, indices_of_points_to_be_used);
	return  model_selection::aic(nll, k);

}

fp_type model_selection::mdl(const commonutil::DataSet& dataset, const Parameters& gmmdesc, const std::vector< idx_type > & indices_of_points_to_be_used)
{
	unsigned int n = indices_of_points_to_be_used.empty() ? dataset.points.cols() : indices_of_points_to_be_used.size();
	unsigned int k = gmmdesc.weights.size();
	unsigned int d = dataset.points.rows();
	fp_type nll = gmmutil::nll(dataset, gmmdesc, indices_of_points_to_be_used);
	return model_selection::mdl(nll, n, d, k);
}

