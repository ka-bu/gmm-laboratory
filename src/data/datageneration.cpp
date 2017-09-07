#include "datageneration.h"

#include "../base.h"
#include "../base/parameters.h"
#include "../base/linalgutil.h"
#include "../base/commonutil.h"
#include "../algorithm/gmm/utils/gmmutil.h"
#include "../settings/settings.h"
#include "utils/gaussian.h"
#include "utils/gaussianmixture.h"
#include "utils/datautil.h"
#include "../ioutil/data_ioutil.h"

#include <stdio.h>



Parameters datageneration::generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						   fp_type separation,
						   fp_type minSqrtEW, fp_type maxSqrtEW, 
						   fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						   fp_type weightExp,
						   fp_type minSqrtEWExp)
{
	
	if(minSqrtEWExp < 0)
	  gmmlab_throw("Use minSqrtEWExp > 0.");
  
	Parameters truth = datageneration::generateRandomGMMWithUniformMeans(d,k,gen,0,maxSqrtEW*maxSqrtEWProportion,minSqrtEW,maxSqrtEW,minSqrtEWProportion,maxSqrtEWProportion,weightExp,minSqrtEWExp);
	
	fp_type separation_tmp  = gmmutil::getSeparation(truth);
	truth.means = separation/separation_tmp*truth.means;
	
// 	if (commonSettings().verbose)
// 		std::cout << "generated " << truth << std::endl;
		
	return truth;
}


Parameters datageneration::generateRandomGMMWithUniformMeans(idx_type d, idx_type k, std::mt19937& gen,
						   fp_type minMean, fp_type maxMean,
						   fp_type minSqrtEW, fp_type maxSqrtEW, 
						   fp_type minSqrtEWProportion, fp_type maxSqrtEWProportion,
						   fp_type weightExp,
						   fp_type minSqrtEWExp)
{
	
	if(minSqrtEWExp < 0)
	  gmmlab_throw("Use minSqrtEWExp > 0.");
  
	Parameters truth;
	
	std::uniform_real_distribution<> urd(0, 1);	

	// (1) create means uniformly at random in a given area
	Vector mean = Vector(d);
	truth.means = Matrix::Zero(d,k);
	for(idx_type l=0; l<k; ++l)
	{	
		commonutil::fill(mean, gen, minMean, maxMean);
		truth.means.col(l) = mean;
	}
		
	// (2) create random covariances with restricted eigenvalues
		
	// create factors that determine the growth of the minimum eigenvalues
	Vector minSqrtEWExpVector = datageneration::getExpIncrWeights(k, minSqrtEWExp);
	minSqrtEWExpVector = minSqrtEWExpVector * 1./minSqrtEWExpVector.maxCoeff();
	// shuffle these factors (unless all factors = 1/k) to avoid a 
	// correlation between minimum eigenvalues and weights (which are determined
	// deterministically
	if(minSqrtEWExp != 0)
	{
		std::vector<fp_type> shuffleVector(minSqrtEWExpVector.size());
		for(idx_type i=0; i<shuffleVector.size(); ++i)
			shuffleVector.at(i) = minSqrtEWExpVector(i);
		std::shuffle<>(shuffleVector.begin(), shuffleVector.end(), gen);
		for(idx_type i=0; i<shuffleVector.size(); ++i)
			minSqrtEWExpVector(i) = shuffleVector.at(i);
	}
	// create covariance matrizes
	Vector eigenvalues = Vector(d);
	Vector traces = Vector(k);
	Matrix randMatrix = Matrix(d,d);
	Matrix randOrthonormalMatrix = Matrix(d,d);
	for(idx_type l=0; l<k; ++l)
	{			
		// create random matrix M
		for (idx_type i=0; i<d; ++i)
			for (idx_type j=0; j<d; ++j)
				randMatrix(i,j)=urd(gen);
			
		// obtain random orthonormal matrix Q by the QR-decomposition of the random matrix M
		Eigen::HouseholderQR<Matrix> qr(randMatrix);
		randOrthonormalMatrix = qr.householderQ();
		
		// create random eigenvalues
		fp_type minEW, maxEW;
		minEW = minSqrtEW + minSqrtEWExpVector(l) * (maxSqrtEW - minSqrtEW);	
		maxEW = (minSqrtEWProportion + urd(gen) * (maxSqrtEWProportion-minSqrtEWProportion))*minEW;
		minEW = pow(minEW, 2);
		maxEW = pow(maxEW, 2);
		commonutil::fill(eigenvalues, gen, minEW, maxEW);
		// make sure that minEW and maxEW appear in the EW's
		eigenvalues(0) = minEW;
		eigenvalues(d-1) = maxEW;
		// Note that it is not necessary to shuffle the eigenvalues, since we multiply the resulting
		// diagonal matrix containing these eigenvalues with a completely random matrix (Q).
		
		traces(l) = eigenvalues.sum();
		
		// random covariance = Q^T * Diag(eigenvalues) * Q
		truth.covariances.push_back(randOrthonormalMatrix.transpose()*eigenvalues.asDiagonal()*randOrthonormalMatrix);
	}
	
	// (3) weights
	truth.weights = datageneration::getExpIncrWeights(k, weightExp);
	// (4) kappas
	truth.kappa = Vector::Constant(k, 1);
	// print some infos	
// 	std::cout << "trace/weight = " << (truth.weights.cwiseQuotient(traces)).transpose() << std::endl;
// 	std::cout << "D= " << truth.means.rows() << std::endl;
	
	return truth;
}

Parameters datageneration::generateRandomGMM(idx_type d, idx_type k, std::mt19937& gen, unsigned int weightExp, fp_type spread)
{
	Vector metaMean = Vector::Zero(d);
	if (commonSettings().verbose)
		std::cout << "the meta mean is " << metaMean.transpose() << std::endl;

	std::normal_distribution<> nd(0, spread*sqrt(fp_type(k)/fp_type(d)));
	Matrix m(d,d);
	for (idx_type i=0; i<d; ++i)
		for (idx_type j=0; j<d; ++j)
			m(i,j)=nd(gen);
	Matrix metaCovariance = m*m.transpose();
	if (commonSettings().verbose)
		std::cout << "the meta covariance is" << std::endl << metaCovariance << std::endl;

	//SingleGaussian metaGauss(metaMean, metaCovariance);
	std::uniform_real_distribution<> urd(0, 1);
	if (k!=d || spread>1)
		nd = std::normal_distribution<>(0, 1);
	Parameters truth;

	std::vector<fp_type> w(k);
	fp_type sum = 0.0;
	for (size_t i=0; i<k; ++i)
	{
		w[i] = pow(urd(gen), weightExp);
		sum += w[i];
	}

	truth.weights = Vector(k);
	truth.means = Matrix(d,k);
	for (size_t i=0; i<k; ++i)
	{
		truth.weights[i] = w[i]/sum;
		//truth.means.col(i) = metaGauss.draw(gen);
		truth.means.col(i) = datautil::drawFromSG(metaMean, metaCovariance, gen, 1);
		for (idx_type i=0; i<d; ++i)
			for (idx_type j=0; j<d; ++j)
				m(i,j)=nd(gen);
		truth.covariances.push_back(m*m.transpose());
		truth.kappa[i] = 1;
	}

	if (commonSettings().verbose)
		std::cout << "generated " << truth << std::endl;

	return truth;
}


void datageneration::generateInputFromGMM(commonutil::DataSet& target, Parameters const& desc, idx_type n, std::mt19937& gen)
{
	assert(desc.means.cols()==desc.weights.size() && desc.covariances.size()==desc.weights.size());

	const idx_type d = desc.means.rows();
	if (desc.means.size()==0)
		gmmlab_throw("Cannot draw data from an empty gmmdesc.");

	if(commonSettings().verbose)
	{
		std::cout << "\n\ngenerating data from gmm..."<< std::endl;
		//std::cout << "  weights = " << (truth.weights).transpose() << std::endl;
		std::cout << "  points per cluster = ";
		for(idx_type i=0; i<desc.components(); ++i)
			printf("%u, ",(unsigned int) std::round(desc.weights(i)*n));
// std::cout << "  parameters = " << desc << std::endl;
		std::cout << std::endl << std::endl;
	}
	
	// resize the input matrix so that each data point will form one column
	target.points.resize(d, n);
	target.weights.resize(n);

	GaussianMixture gm(desc);
	for (std::size_t i=0; i<n; ++i)
	{
		target.weights[i] = 1;
		
		target.points.col(i) = gm.draw(gen);
		fp_type density = gm.density(target.points.col(i));
		if (density<=0)
			std::cout << "datagenutil::generateInputFromGMM - WARNING!!! Drew point " << target.points.col(i) << " having density " << density << "." << std::endl;
			
	}
	
	
DataIOUtil::printStatistis(target);
	
}





void datageneration::generateCorrelatedNoisyInputFromGMM(commonutil::DataSet& target, Parameters const& desc, idx_type n, fp_type sqrtEWFraction, std::mt19937& gen)
{
	assert(desc.means.cols()==desc.weights.size() && desc.covariances.size()==desc.weights.size());

	const idx_type d = desc.means.rows();
	const idx_type k = desc.means.cols();
	
	if(sqrtEWFraction <= 0)
		gmmlab_throw(" SqrtEWFraction should be > 0.");
	if (d==0 || k==0)
		gmmlab_throw(" Cannot draw data from an empty gmmdesc.");

	// resize
	target.points.resize(d, n);
	target.weights.resize(n);

	// determine noise per component (centered around zero!)
	std::vector<SingleGaussian> noisePerComponent; 
	Vector zero = Vector::Zero(d);
	for(idx_type i=0; i<k; ++i)
	{
		Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(desc.covariances.at(i));
		if (eigensolver.info() != Eigen::Success)
			gmmlab_throw(" Eigensolver failed.");
		Vector eigenvalues = eigensolver.eigenvalues();
		eigenvalues *= pow(sqrtEWFraction,2);		
		Matrix q = eigensolver.eigenvectors();
		SingleGaussian gaussian(zero, q*eigenvalues.asDiagonal()*q.transpose());
		noisePerComponent.push_back(gaussian);
	}
	
	// draw data and "add" gaussian error
	GaussianMixture gm(desc);
	for (idx_type i=0; i<n; ++i)
	{
		// draw a point according to the truth
		commonutil::CompleteData cd = gm.drawCompleteData(gen);
		fp_type density = gm.density(cd.point);
		if (density<=0)
			std::cout << "datagenutil::generateInputFromGMM - WARNING!!! Drew point " << target.points.col(i) << " having density " << density << "." << std::endl;
		
		// draw a point according to noise distribution, which is centered around zero, then shift it to the 
		// point drawn according to the truth
		Vector noisyPoint = noisePerComponent.at(cd.source).draw(gen);
		noisyPoint += cd.point;
		
		target.points.col(i) = noisyPoint;
		target.weights[i] = 1;
	}
}


void datageneration::generateDifficultNoisyInputFromGMM(commonutil::DataSet& target, Parameters const& desc, idx_type n, std::mt19937& gen)
{
	assert(desc.means.cols()==desc.weights.size() && desc.covariances.size()==desc.weights.size());

	const idx_type d = desc.means.rows();
	const idx_type k = desc.means.cols();
	
	if (d==0 || k==0)
		gmmlab_throw(" Cannot draw data from an empty gmmdesc.");

	// resize
	target.points.resize(d, n);
	target.weights.resize(n);

	// determine minEw/maxEW over *all* covariances
	fp_type minEW;
	fp_type maxEW;
	for(idx_type i=0; i<k; ++i)
	{
		Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(desc.covariances.at(i));
		if (eigensolver.info() != Eigen::Success)
			gmmlab_throw(" Eigensolver failed.");
		Vector eigenvalues = eigensolver.eigenvalues();
		if(i==0)
		{
			minEW = eigenvalues.minCoeff();
			maxEW = eigenvalues.maxCoeff();
		}
		else
		{
			minEW = std::min(minEW, eigenvalues.minCoeff());
			maxEW = std::max(maxEW, eigenvalues.maxCoeff());
		}
	}
	
	// noisePerComponent (dependent on the *overall* minimum/maximum eigenvalue)
	std::vector<SingleGaussian> noises;
	Vector densities = datageneration::getExpIncrWeights(k, 1.);
	Vector zero = Vector::Zero(d);
	for(idx_type i=k; i>0; --i)
	{
		SingleGaussian gaussian(zero, (minEW + (1.*i)/(1.*k) * ( maxEW - minEW ))*Matrix::Identity(d,d) );
		noises.push_back(gaussian);
	}
	
	
	// draw data and "add" gaussian error
	GaussianMixture gm(desc);
	for (idx_type i=0; i<n; ++i)
	{
		// draw a point according to the truth
		commonutil::CompleteData cd = gm.drawCompleteData(gen);
		fp_type density = gm.density(cd.point);
		if (density<=0)
			std::cout << "datagenutil::generateInputFromGMM - WARNING!!! Drew point " << target.points.col(i) << " having density " << density << "." << std::endl;
		
		// draw a point according to noise distribution, which is centered around zero, then shift it to the 
		// point drawn according to the truth
		idx_type noise = commonutil::randomIndex<>(densities, gen);
		Vector noisyPoint = noises.at(noise).draw(gen);
		noisyPoint += cd.point;
		
		target.points.col(i) = noisyPoint;
		target.weights[i] = 1;
	}
	
}

void datageneration::addUniformNoise(commonutil::DataSet& dataset, idx_type n, std::mt19937& gen)
{	
	const idx_type old_n = dataset.points.cols();
	const idx_type d = dataset.points.rows();
	
	if(old_n==0||d==0)
		gmmlab_throw("Given dataset must not be empty.");
	
	const BoundingBox boundingBox = datageneration::boundingBox(dataset);
	Vector anchor = boundingBox.min;
	Vector interval = boundingBox.max-boundingBox.min;
	anchor = anchor - 0.1*interval;
	interval = interval + 0.2*interval;
	
	datageneration::addUniformBox(dataset, d, n, gen, anchor, interval);
}


void datageneration::addUniformProjectedToSubspaceNoise(commonutil::DataSet& dataset, idx_type n, std::mt19937& gen, bool verbose)
{
	const idx_type old_n = dataset.points.cols();
	const idx_type d = dataset.points.rows();
	
	if(old_n==0||d==0)
		gmmlab_throw("Given dataset must not be empty.");
		
	// compute orthonormal basis of a random subspace with >=1 and <=d dimensions
	Matrix m = Matrix(d,d);
	commonutil::fill(m, gen, 0, 1);
	Eigen::HouseholderQR<Matrix> qr(m);
	Matrix q = qr.householderQ();
	std::uniform_int_distribution<> uid(1,d-1);
	unsigned int subd = uid(gen);
	q.conservativeResize(d, subd);
	q = q*q.transpose(); // projection
	
// 	std::cout << " Q*Q^T = " << q << std::endl;
// 	std::cout << "subd = " << subd << std::endl;
	
	// resize dataset
	dataset.points.conservativeResize(d, old_n+n);
	dataset.weights.conservativeResize(old_n+n);

	// draw from the bounding box and project it to the subspace
	const BoundingBox boundingBox = datageneration::boundingBox(dataset);
	
// 	std::cout << " BoundingBox = " << boundingBox.min.transpose() << " and " << boundingBox.max.transpose() << std::endl;
	
	Vector tmp = Vector(d);
	for(idx_type i=0; i<n; ++i)
	{
		commonutil::fill(tmp, gen, boundingBox.min, boundingBox.max);
		
// 		std::cout << "tmp = " << tmp.transpose() << std::endl;
// 		std::cout << "tmp = " << (q*tmp).transpose() << std::endl;
		
		dataset.points.col(old_n+i) = q*tmp;
		dataset.weights(old_n+i) = 1.;		
	}
	
	if(verbose)
		std::cout << "datagenutil::addUniformProjectedToSubspaceNoise() - Added noise in " << subd << " dimensions." << std::endl;
}



void datageneration::addUniformShapes(commonutil::DataSet& dataset, idx_type d, idx_type n, std::mt19937& gen, idx_type k, fp_type weightFloatingExp,
							fp_type minMean, fp_type maxMean, fp_type minDiameter, fp_type maxDiameter)
{	
	if(n<=0)
		return;
	
	unsigned int differentShapes = 4;
	std::uniform_int_distribution<> uid (0,differentShapes-1);
	std::uniform_real_distribution<> urd(0,1);
	
	Vector tmp1 = Vector(d);
	Vector tmp2 = Vector(d);
	
	//Vector pointsPerShape = Vector::Constant(shapes, (unsigned int) round(((unsigned int) n)/(1.*shapes)));
	Vector pointsPerShape = datageneration::getExpIncrNumbers(k, weightFloatingExp, n);
	
	fp_type minRadius, maxRadius;
	fp_type minArcThickness = (maxDiameter-minDiameter)/4.;
	
	for(idx_type i=0; i<k; ++i)
	{
		unsigned int shape = uid(gen);
		
		switch(shape){
			case 0: // Bridge
				commonutil::fill(tmp1, gen, minMean, maxMean); // anchor
				commonutil::fill(tmp2, gen, minDiameter, maxDiameter); // interval
				datageneration::addUniformBridge(dataset, d,pointsPerShape(i), gen, tmp1, tmp2);
				break;
			case 1: // Arc
				commonutil::fill(tmp1, gen, minMean, maxMean); // center
				minRadius = minDiameter/2.+(maxDiameter-minDiameter)/2.*urd(gen);
				maxRadius = minRadius+(maxDiameter/2.-minRadius)*urd(gen);
				if(maxRadius-minRadius < minArcThickness)
					maxRadius = minRadius + minArcThickness;
				datageneration::addUniformArc(dataset, d, pointsPerShape(i), gen, tmp1, minRadius, maxRadius);
				break;
			case 2: // Box
				commonutil::fill(tmp1, gen, minMean, maxMean); // anchor
				commonutil::fill(tmp2, gen, minDiameter, maxDiameter); // interval
				datageneration::addUniformBox(dataset, d, pointsPerShape(i), gen, tmp1, tmp2);
				break;
			default: // Sphere
				commonutil::fill(tmp1, gen, minMean, maxMean); // center
				minRadius = minDiameter/2.+(maxDiameter-minDiameter)/2.*urd(gen);
				maxRadius = minRadius+(maxDiameter/2.-minRadius)*urd(gen);
				datageneration::addUniformSphere(dataset, d, pointsPerShape(i), gen, tmp1, minRadius, maxRadius);
				break;
		}
	}
}



void datageneration::addUniformArc(commonutil::DataSet& dataset, idx_type d,  idx_type n, std::mt19937& gen, Vector center, fp_type minRadius, fp_type maxRadius)
{	
	const idx_type old_n = dataset.points.cols();
	
	if(old_n!=0 && dataset.points.rows()!=d)
		gmmlab_throw("d should be equal to the dimension of the points in the given data set.");
	if(minRadius<0 || minRadius>=maxRadius)// || minFraction < -1 || minFraction >= maxFraction || maxFraction>1)
		gmmlab_throw("Invalid radii.");
	if(center.size() != d)
		gmmlab_throw("Invalid center.");
	if(d<2)
		gmmlab_throw("Dimension should be at least 2");
		
	// resize
	dataset.points.conservativeResize(d,old_n+n);
	dataset.weights.conservativeResize(old_n+n);
	
	// pick 2 dimensions uniformly at random
	// restricted to these dimensions, the points will form an arc
	std::uniform_int_distribution<> uid(0,d-1);
	Vector arcDimensions = Vector(2);
	arcDimensions << uid(gen), uid(gen);

	std::uniform_real_distribution<> urd(0,1);
	for(idx_type i=0; i<n; ++i)
	{
		Vector point = Vector(d);
		// sample entries in the arcDimensions
		Vector arcPoint = Vector(2);
		fp_type norm = 0;
		while(norm==0)
		{
			arcPoint << (urd(gen)), (2*urd(gen)-1);
			norm = arcPoint.norm();
		}
		arcPoint = (minRadius + (maxRadius-minRadius)*urd(gen))* arcPoint/norm;
		// sample entries for all dimensions and overwrite those entries in the arcDimensions
		commonutil::fill(point, gen, minRadius, maxRadius);
		for(idx_type j=0; j<2; ++j)
			point(arcDimensions(j)) = arcPoint(j);
		
		dataset.points.col(n+i) = point+center;
		dataset.weights(n+i) = 1.;
	}
}

void datageneration::addUniformBridge(commonutil::DataSet& dataset, idx_type d,  idx_type n, std::mt19937& gen, Vector anchor, Vector const& interval)
{
	if(d==1)
		gmmlab_throw("Dimension should be > 2.");
	std::uniform_int_distribution<> uid(0,d-1);
	idx_type d1 = uid(gen);
	idx_type d2 = uid(gen);
	while(d1==d2)
		d2 = uid(gen);
	Vector box = interval;
	box(d1) = box(d1)/3.;
	box(d2) = box(d2)/2.;
	
	unsigned int pointsPerBox =  (unsigned int) round(((unsigned int) n)/6.);
	
	datageneration::addUniformBox(dataset, d, pointsPerBox, gen,  anchor, box);
	anchor(d2) = anchor(d2) + interval(d2)/2.;
	datageneration::addUniformBox(dataset, d,  pointsPerBox, gen, anchor, box);
	anchor(d1) = anchor(d1) + interval(d1)/4.;
	datageneration::addUniformBox(dataset, d, pointsPerBox, gen, anchor, box);
	anchor(d1) = anchor(d1) + interval(d1)/4.;
	datageneration::addUniformBox(dataset, d,  pointsPerBox,gen,  anchor, box);
	anchor(d1) = anchor(d1) + interval(d1)/4.;
	datageneration::addUniformBox(dataset, d,  pointsPerBox,gen,  anchor, box);
	anchor(d2) = anchor(d2) - interval(d2)/2.;
	datageneration::addUniformBox(dataset, d,  pointsPerBox,gen,  anchor, box);
}

void datageneration::addUniformBox(commonutil::DataSet& dataset, idx_type d,  idx_type n,  std::mt19937& gen, Vector const& anchor, Vector const& interval)
{		
	const idx_type old_n = dataset.points.cols();
	
	if(interval.rows()!=d || anchor.rows()!=d || interval.rows()!=d)
		gmmlab_throw("Interval or anchor have the wrong dimension.");
	
	// resize
	dataset.points.conservativeResize(d,old_n+n);
	dataset.weights.conservativeResize(old_n+n);
	
	// add points that are drawn uniformly at random
	std::uniform_real_distribution<> urd(0, 1);	
	for(idx_type i=0; i<n; ++i)
	{
		Vector newPoint = Vector(d);
		for(idx_type j=0; j<d; ++j)
			newPoint(j) = urd(gen)*interval(j);
		dataset.points.col(old_n+i) = newPoint+anchor;
		dataset.weights(old_n+i)=1.;
	}		
}


void datageneration::addUniformSphere(commonutil::DataSet& dataset, idx_type d,  idx_type n, std::mt19937& gen, Vector center, fp_type minRadius, fp_type maxRadius)
{
	const idx_type old_n = dataset.points.cols();
	
	if(n!=0 && dataset.points.rows()!=d)
		gmmlab_throw(" d should be equal to the dimension of the points in the given data set.");
	if(minRadius<0 || minRadius>=maxRadius)// || minFraction < -1 || minFraction >= maxFraction || maxFraction>1)
		gmmlab_throw(" Invalid radii.");
	if( center.size() != d)
		gmmlab_throw(" Invalid center.");
	
	// resize
	dataset.points.conservativeResize(d,old_n+n);
	dataset.weights.conservativeResize(old_n+n);
		
	std::uniform_real_distribution<> urd(0, 1);
	for(idx_type i=0; i<n; ++i)
	{		

		Vector point = Vector(d);
		fp_type norm = 0;
		while(norm==0)
		{
			commonutil::fill(point, gen, -1,1);
			norm = point.norm();
		}
		point = (minRadius + (maxRadius-minRadius)*urd(gen))* point/norm;
		// shift the points by center
		dataset.points.col(old_n+i) = point+center;
		dataset.weights(old_n+i) = 1;
	}
	
			// uniformly choose a radius in [minRadius, maxRadius], 
		// then pick a point from the resulting sphere (with mean 0) 
		// uniformly at random
	// 	// following does *not* distribute points uniformly on the sphere...
	// 	for(idx_type p=0; p<number; ++p)
	// 	{
	// 		fp_type r = minRadius + (maxRadius-minRadius)*urd(gen);
	// 		Vector phi = Vector(d);
	// 		Vector point = Vector(d);
	// 		fp_type diffAngle = maxAngle-minAngle;	
	// 		for(idx_type i=0; i<d; ++i)
	// 		{
	// 			if(i<d-1)
	// 			{
	// 				if(i==d-2)
	// 					phi(i) = urd(gen)*2*M_PI;
	// 				else
	// 					phi(i) = minAngle + diffAngle*urd(gen);
	// 				point(i) = r*cos(phi(i));
	// 			}
	// 			else
	// 				point(d-1) = r*sin(phi(d-2));
	// 			for(idx_type j=0; j<i; ++j)
	// 				point(i) *= sin(phi(j));
	// 		}
	// 				
	// 		dataset.points.col(n+p) = point;
	// 		dataset.weights(n+p) = 1;
	// 	}
}
	

Vector datageneration::getExpIncrNumbers(idx_type k, fp_type exponent, idx_type number)
{
	Vector output = datageneration::getExpIncrWeights(k, exponent);
	output = output*number;
	for(idx_type i=0; i<k-1; ++i)
		output(i) = round(output(i));
	output(k) = 0;
	output(k) = number - output.sum();
	return output;
}

Vector datageneration::getExpIncrWeights(idx_type k, fp_type exponent)
{
	Vector output = Vector(k);
	fp_type factor = pow(2,exponent);
	output(0) = 1.;
	for(idx_type i=1; i<k; ++i)
	  output(i) = output(i-1)*factor;	
	output = output/output.sum();
	return output;
}

datageneration::BoundingBox datageneration::boundingBox(commonutil::DataSet const& dataset)
{
	const idx_type n = dataset.points.cols();
	const idx_type d = dataset.points.rows();
	
	if(n==0||d==0)
		gmmlab_throw("Given dataset must not be empty.");
		
	BoundingBox boundingBox;
	boundingBox.min = dataset.points.col(0);
	boundingBox.max = boundingBox.min;
	for(idx_type i=1; i<n; ++i)
	{
		boundingBox.min = boundingBox.min.cwiseMin(dataset.points.col(i));
		boundingBox.max = boundingBox.max.cwiseMax(dataset.points.col(i));
	}
	return boundingBox;
}
