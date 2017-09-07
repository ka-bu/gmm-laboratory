#include "commonutil.h"


void commonutil::get_uniform_subset(commonutil::DataSet const& dataset, commonutil::DataSet& samples, idx_type number, std::mt19937& gen)
{

	idx_type n = dataset.points.cols();
	idx_type d = dataset.points.rows();
		
	std::uniform_int_distribution<> uid(0,n-1);
	
	samples.points.resize(d,number);
	samples.weights.resize(number);
	
	for(idx_type j=0; j<number; ++j)
	{
		idx_type index = uid(gen);		
		samples.points.col(j) = dataset.points.col(index);
		samples.weights(j) = dataset.weights(index);
	}
}


std::vector< idx_type > commonutil::sort_indices(const Vector& vector, bool increasing)
{
	idx_type n = vector.size();
	std::vector<idx_type> idx(n);
	for(idx_type i=0; i<n; ++i)
		idx[i] = i;
	std::sort(idx.begin(),idx.end(),
		  [&vector, increasing](idx_type i1, idx_type i2){
			  if(increasing)
				  return vector(i1) < vector(i2);
			  else
				  return vector(i1) > vector(i2);
		});
	return idx;
}

std::vector< idx_type > commonutil::nth_element_indices(const Vector& vector, idx_type n, bool increasing)
{
	idx_type vector_size = vector.size();
	std::vector<idx_type> idx(vector_size);
	for(idx_type i=0; i<vector_size; ++i)
		idx[i] = i;
	std::nth_element(idx.begin(),idx.begin()+n-1,idx.end(),
		  [&vector, increasing](idx_type i1, idx_type i2){
			  if(increasing)
				  return vector(i1) < vector(i2);
			  else
				  return vector(i1) > vector(i2);
		});
	return idx;
}


void commonutil::erase(Vector& v, idx_type i)
{
	const idx_type s = v.size();
	if(s==0)
		gmmlab_throw("The vector is empty, there is no entry to erase.");
	if(i<0 || i>=s)
		gmmlab_throw("There does not exist an entry with this index.");
	if(i < s-1)
		std::copy(v.data() + i + 1, v.data() + v.size(), v.data() + i);
    v.conservativeResize(s-1);
}


void commonutil::eraseColumn(Matrix& m, idx_type i)
{
	const idx_type r = m.rows();
	const idx_type c = m.cols();
	if(c==0)
		gmmlab_throw("The matrix is empty, there is no column to erase.");
	if(i<0 || i>=c)
		gmmlab_throw("There does not exist a column with this index.");
	if(i<c-1)
		// shift columns to the left. note: (c-1) is the number of the last column
		m.block(0,i,r,(c-1)-(i+1)+1) = m.block(0,i+1,r,(c-1)-(i+1)+1);
	m.conservativeResize(Eigen::NoChange, c-1);
}


bool commonutil::componentWiseDifferent(Vector p1, Vector p2)
{
	if(p1.size()!=p2.size())
		gmmlab_throw("Vectors have different size");
	for(idx_type i=0; i<p1.size(); ++i){
		if(p1(i)==p2(i))
			return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& os, commonutil::DataSet const& data)
{
	const idx_type s = data.weights.size();

	assert(data.points.cols()==s);

	os << "data set with " << s << " points" << std::endl;
	os << data.points << std::endl << std::endl
		<< data.weights.transpose() << std::endl << std::endl;
	return os;
}

void commonutil::appendData(commonutil::DataSet& dataset, commonutil::DataSet const& additionalData)
{
	const idx_type n = dataset.points.cols();
	const idx_type d = dataset.points.rows();
	const idx_type add_n = additionalData.points.cols();
	const idx_type add_d = additionalData.points.rows();
	
	if(d != add_d)
		gmmlab_throw("")
	
	dataset.points.conservativeResize(d,n+add_n);
	dataset.weights.conservativeResize(n+add_n);
}


void commonutil::fill(Matrix& matrix, std::mt19937& gen, fp_type min, fp_type max)
{
	const idx_type d = matrix.rows();
	const idx_type n = matrix.cols();
	if(d<=0 || n<=0)
		gmmlab_throw("Matrix should have dimensions > 0.")
	std::uniform_real_distribution<>urd(min,max);
	for(idx_type i=0; i<d; ++i)
		for(idx_type j=0; j<n; ++j)
			matrix(i,j) = urd(gen);
}
	
void commonutil::fill(Vector& vector, std::mt19937& gen, fp_type min, fp_type max)
{
	const idx_type d = vector.size();
	if(d<=0)
		gmmlab_throw("Vector should have dimension > 0.")
	std::uniform_real_distribution<>urd(min,max);
	for(idx_type i=0; i<d; ++i)
		vector(i) = urd(gen);
}

void commonutil::fill(Vector& vector, std::mt19937& gen, Vector min, Vector max)
{
	const idx_type d = vector.size();
	if(d<=0)
		gmmlab_throw("Vector should have dimension > 0.")
	std::uniform_real_distribution<>urd(0,1);
	for(idx_type i=0; i<d; ++i)
		vector(i) = min(i) + urd(gen)*(max(i)-min(i));
}

void commonutil::map_indices ( std::vector<idx_type>& indices, const std::vector<idx_type>& mapping )
{
	if(*std::max_element(indices.begin(), indices.end())>= mapping.size())
		gmmlab_throw("Mapping does not suffice.");
	for(idx_type i=0; i<indices.size(); ++i)
		indices.at(i) = mapping.at(indices.at(i));
}


