#include "../base.h"

#include "gmmdisplay.h"

#include <iostream>
#include <cstdio>
#include <AntTweakBar.h>
#include <ctime>
#include <boost/math/distributions/chi_squared.hpp>

#ifndef GLFW_INCLUDE_GLU
#define GLFW_INCLUDE_GLU
#endif
#include <GLFW/glfw3.h>

#include "../gmmlab.h"
#include "../base/linalgutil.h"
#include "../ioutil/gmm_ioutil.h"
#include "../algorithm/kmeans/utils/kmeansutil.h"
#include "../algorithm/gmm/utils/gmmutil.h"
#include "../data/datageneration.h"


GMMDisplay::GMMDisplay(Parameters& t, commonutil::DataSet const& ip,
	std::vector<commonutil::DataSet> const& ss, std::vector<Parameters> const& slt,
	std::vector<std::string> const& tg, std::vector<double> const& rt,
	std::vector<fp_type> const& nll)
	: truth(t), input(ip), sampleSets(ss), solutions(slt), tags(tg), runtimes(rt), nlls(nll), csd(ip.points.rows())
{
	this->statistics.truthNLL = 0;
	this->statistics.separation = -1;
	this->statistics.tag = "NONE";
	this->statistics.nll = 0;
	this->statistics.runtime = 0;
	this->alpha = 0;
	this->quantile = 1;
	this->pointSize = displaySettings().pointSize;
	
	initGLFW();
	initGL();

	prepareDisplayLists();

	idx_type num = this->input.points.rows();
	observationBar->setProjectionRange(1,num);
	
	std:size_t size = this->solutions.size();
	observationBar->setSolutionRange(0,size-1);
	if (size>0)
		observationBar->setComponentRange(0,this->solutions[0].weights.size());

	// compute cost of truth
	if (!t.empty())
	{
	  if(commonSettings().costmeasure == GMMLab::COST_GAUSSIAN)
	    this->statistics.truthNLL = gmmutil::nll(ip,truth);
	  else if(commonSettings().costmeasure == GMMLab::COST_KMEANS)
	    this->statistics.truthNLL = kmeansutil::kmeanscost(ip, truth.means);
	  else 
		gmmlab_throw("unkown cost function")
	  
	  this->statistics.separation = gmmutil::getSeparation(truth);
	}

	resetProjection();
}

GMMDisplay::~GMMDisplay()
{
	// Terminate AntTweakBar and GLFW
	TwTerminate();
	glfwTerminate();
}

// Callback function called by GLFW when window size changes
void GMMDisplay::WindowSizeCB(GLFWwindow* window, int width, int height)
{
    // Set OpenGL viewport and camera
    glViewport(0, 0, width, height);
    displaySettings().aspectRatio = (double)width/height;

    // Send the new window size to AntTweakBar
    TwWindowSize(width, height);

//     std::cout << "new window size: width=" << width << " height=" << height << std::endl;
//     std::cout << "aspect ratio = " << displaySettings().aspectRatio << std::endl;

}

void GMMDisplay::prepareDisplayLists()
{
	this->pointsDL = 0;

	this->hiSphereDL = glGenLists(1);
	glNewList(this->hiSphereDL, GL_COMPILE);
	gluSphere(gluNewQuadric(), 1, 50, 50);
	glEndList();

	this->loSphereDL = glGenLists(1);
	glNewList(this->loSphereDL, GL_COMPILE);
	gluSphere(gluNewQuadric(), 1, 8, 8);
	glEndList();

	this->boxDL = glGenLists(1);
	glNewList(this->boxDL, GL_COMPILE);
	glBegin(GL_QUADS);
	glNormal3f(0,0,-1); glVertex3d(-1,-1,-1); glVertex3d(-1,+1,-1); glVertex3d(+1,+1,-1); glVertex3d(+1,-1,-1); // front face
	glNormal3f(0,0,+1); glVertex3d(-1,-1,+1); glVertex3d(+1,-1,+1); glVertex3d(+1,+1,+1); glVertex3d(-1,+1,+1); // back face
	glNormal3f(-1,0,0); glVertex3d(-1,-1,-1); glVertex3d(-1,-1,+1); glVertex3d(-1,+1,+1); glVertex3d(-1,+1,-1); // left face
	glNormal3f(+1,0,0); glVertex3d(+1,-1,-1); glVertex3d(+1,+1,-1); glVertex3d(+1,+1,+1); glVertex3d(+1,-1,+1); // right face
	glNormal3f(0,-1,0); glVertex3d(-1,-1,-1); glVertex3d(+1,-1,-1); glVertex3d(+1,-1,+1); glVertex3d(-1,-1,+1); // bottom face
	glNormal3f(0,+1,0); glVertex3d(-1,+1,-1); glVertex3d(-1,+1,+1); glVertex3d(+1,+1,+1); glVertex3d(+1,+1,-1); // top face
	glEnd();

	glBegin(GL_LINES);
	glVertex3d(-.3,0,0); glVertex3d(.3,0,0);
	glVertex3d(0,-.3,0); glVertex3d(0,.3,0);
	glVertex3d(0,0,-.3); glVertex3d(0,0,.3);
	glEnd();
	glEndList();
}

GLdouble* GMMDisplay::convertMatrix(Matrix const& in)
{
	Matrix m = Matrix::Zero(3,3);
	idx_type r = in.rows();
	if (r>3)
		r = 3;
	idx_type c = in.cols();
	if (c>3)
		c = 3;
	m.block(0,0,r,c) = in.block(0,0,r,c);
	
	GLdouble* out = new GLdouble[16];

	out[0*4+0] = GLdouble(m(0,0));
	out[1*4+0] = GLdouble(m(0,1));
	out[2*4+0] = GLdouble(m(0,2));

	out[0*4+1] = GLdouble(m(1,0));
	out[1*4+1] = GLdouble(m(1,1));
	out[2*4+1] = GLdouble(m(1,2));

	out[0*4+2] = GLdouble(m(2,0));
	out[1*4+2] = GLdouble(m(2,1));
	out[2*4+2] = GLdouble(m(2,2));

	out[0*4+3] = GLdouble(0);
	out[1*4+3] = GLdouble(0);
	out[2*4+3] = GLdouble(0);

	out[3*4+0] = GLdouble(0);
	out[3*4+1] = GLdouble(0);
	out[3*4+2] = GLdouble(0);

	out[15] = GLdouble(1);

	return out;
}

// Routine to convert a quaternion to a 4x4 matrix
// ( input: quat = float[4]  output: mat = float[4*4] )
void GMMDisplay::convertQuaternionToMatrix(const double *quat, double *mat)
{
    double yy2 = 2.0f * quat[1] * quat[1];
    double xy2 = 2.0f * quat[0] * quat[1];
    double xz2 = 2.0f * quat[0] * quat[2];
    double yz2 = 2.0f * quat[1] * quat[2];
    double zz2 = 2.0f * quat[2] * quat[2];
    double wz2 = 2.0f * quat[3] * quat[2];
    double wy2 = 2.0f * quat[3] * quat[1];
    double wx2 = 2.0f * quat[3] * quat[0];
    double xx2 = 2.0f * quat[0] * quat[0];
    mat[0*4+0] = - yy2 - zz2 + 1.0f;
    mat[0*4+1] = xy2 + wz2;
    mat[0*4+2] = xz2 - wy2;
    mat[0*4+3] = 0;
    mat[1*4+0] = xy2 - wz2;
    mat[1*4+1] = - xx2 - zz2 + 1.0f;
    mat[1*4+2] = yz2 + wx2;
    mat[1*4+3] = 0;
    mat[2*4+0] = xz2 + wy2;
    mat[2*4+1] = yz2 - wx2;
    mat[2*4+2] = - xx2 - yy2 + 1.0f;
    mat[2*4+3] = 0;
    mat[3*4+0] = mat[3*4+1] = mat[3*4+2] = 0;
    mat[3*4+3] = 1;
}

void GMMDisplay::drawGaussian(Vector const& mean, Matrix covariance, float const* rgba)
{
	glPushMatrix();
	
	Vector p = this->projection*mean;
	glTranslated(p[0], p[1], p[2]);

	Matrix covchol;

	if(observationSettings().meansOnly > 0)
	{
		idx_type d = covariance.rows();
		covariance = Matrix::Identity(d,d)*observationSettings().meansOnly;
	}
		
	if(covariance.rows() >= 3)
	{
		covchol = this->projection*(quantile*covariance)*this->projection.transpose();
		linalg::cholesky(covchol);
	}
	else if(covariance.rows() == 1)
	{
		covchol = sqrt(covariance(0,0))*Matrix::Identity(3,3);
	}
	else if(covariance.rows() == 2)
	{;
		covchol = Matrix::Identity(3,3);
		covchol.block(0,0,2,2) = quantile*covariance;
		linalg::cholesky(covchol);
	}

	
	GLdouble* m = convertMatrix(covchol);
	glMultMatrixd(m);
	delete[] m;

	glColor4fv(rgba);

	if( displaySettings().wireframe )
	{
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		if (this->loSphereDL!=0)
			glCallList(this->loSphereDL);
		else
			gluSphere(gluNewQuadric(), 1, 8, 8);
	}
	else
	{
		glEnable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		if (this->hiSphereDL!=0)
			glCallList(this->hiSphereDL);
		else
			gluSphere(gluNewQuadric(), 1, 32, 32);
	}

	glPopMatrix();
}

void GMMDisplay::drawGMM(Parameters const& desc, float const* normalColor, float const* selectionColor)
{
	std::size_t k = desc.weights.size();
	fp_type maxWeight = desc.weights.maxCoeff();
	for (std::size_t i=0; i<k; ++i)
	{
//		std::cout << "trying to draw gaussian " << i << std::endl;
//		float mixing = (1+4*desc.weights[i])/5;
		
		float mixing = desc.weights[i]/maxWeight;
		float weightedColor[4];
		weightedColor[0] = normalColor[0]*mixing;
		weightedColor[1] = normalColor[1]*mixing;
		weightedColor[2] = normalColor[2]*mixing;
		weightedColor[3] = normalColor[3];
		
		const float* col;
		if(i+1==observationSettings().component&&selectionColor!=NULL)
			col = selectionColor;
		else if(observationSettings().meansOnly>0)
			col = normalColor;
		else 
			col = weightedColor;
		
		drawGaussian(desc.means.col(i), desc.covariances[i],col);
	}
}

void GMMDisplay::drawPoints(Matrix const& points)
{
	
	std::vector<idx_type> partition;
	if(observationSettings().partition)
	{
		if(commonSettings().costmeasure == GMMLab::COST_KMEANS)
			partition = kmeansutil::kmeansPartition(this->input.points, this->solutions[observationSettings().solution].means);
	}
	
	idx_type k = 0;
	if (!partition.empty())
		k = *std::max_element(partition.begin(), partition.end()) + 1;
	if (k > 24)
	{
		std::cout << "Warning: Can not color more than 24 cluster" << std::endl;
		partition = std::vector<idx_type>();
	}
	idx_type n = points.cols();
	
	Matrix p = this->projection*points;

// 	std::cout << "draw Points... solution = " << observationSettings().solution << std::endl;

	
	glPointSize(displaySettings().pointSize);
	
	if(!partition.empty())
	{
		glDisable(GL_LIGHTING);
	}
	else
		glColor4ubv(displaySettings().pointColor);
	
	glBegin(GL_POINTS);	
	for (idx_type i=0; i<n; ++i)
	{
		if(!partition.empty())
			glColor4fv(displaySettings().clusterColors[partition.at(i)]);
		glVertex3d(p(0,i), p(1,i), p(2,i));
		
	}
	glEnd();

	
	
// 	idx_type i = 15874;
// 	glColor4fv(displaySettings().highlightColor);
// 		glPushMatrix();
// 		glTranslated(p(0,i), p(1,i), p(2,i));
// 		double factor=0.1*pow(2,displaySettings().length);
//  		glScaled(factor,factor,factor);
// 		gluSphere(gluNewQuadric(), 1, 3, 2); 
// 		glPopMatrix();
}

void GMMDisplay::drawBoxes(commonutil::DataSet const& data)
{
	idx_type n = data.points.cols();
	if (n==0)
		return;
	double factor=pow(2,displaySettings().length);

	glColor3fv(displaySettings().boxColor);

	if (displaySettings().wireframe)
	{
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else
	{
		glEnable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	Matrix p = projection*data.points;
	for (size_t i=0; i<n; ++i)
	{
//		fp_type l = factor*pow(data.weights[i],1.0/3.0);
		fp_type l = factor*data.weights[i]; 
		glPushMatrix();
		glTranslated(p(0,i), p(1,i), p(2,i));
		glScaled(l,l,l);

		if (this->boxDL!=0)
			glCallList(this->boxDL);
		else
			gluSphere(gluNewQuadric(), 1, 3, 2); 
			//gluSphere(gluNewQuadric(), 10000., 300, 200); 
		glPopMatrix();
	}
}

void GMMDisplay::adjustDisplay()
{
	std::size_t n = this->input.points.cols();
	if (n == 0)
		return;

	if (this->pointsDL!=0)
		glDeleteLists(this->pointsDL, 1);
	this->pointsDL = glGenLists(1);

	glPointSize(displaySettings().pointSize);
	glNewList(this->pointsDL, GL_COMPILE);	
	drawPoints(this->input.points);
	glEndList();

	Matrix p = projection*this->input.points;

	Vector minEdge = p.rowwise().minCoeff();
	Vector maxEdge = p.rowwise().maxCoeff();
	Vector center = (minEdge+maxEdge)/2;
	this->lightDistance = (maxEdge-minEdge).norm();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	float v[4]; // will be used to set light parameters
	v[0] = v[1] = v[2] = this->lightDistance; v[3] = 0;
	glLightfv(GL_LIGHT0, GL_POSITION, v);

	displaySettings().distance = this->lightDistance;
	displaySettings().zFar = 2.0*displaySettings().distance;
	displaySettings().zNear = displaySettings().zFar / 100.0;

	char defstring[100];
//	std::snprintf(defstring, sizeof(defstring), " display/znear step=%f min=%f ", displaySettings().zNear/10, displaySettings().zNear/10);
	std::sprintf(defstring, " display/znear step=%f min=%f ", displaySettings().zNear/10, displaySettings().zNear/10);
	TwDefine(defstring);
//	std::snprintf(defstring, sizeof(defstring), " display/zfar step=%f min=%f ", displaySettings().zFar/100, displaySettings().zFar/100);
	std::sprintf(defstring, " display/zfar step=%f min=%f ", displaySettings().zFar/100, displaySettings().zFar/100);
	TwDefine(defstring);

	observationSettings().focus[0] = center[0];
	observationSettings().focus[1] = center[1];
	observationSettings().focus[2] = center[2];
}

void GMMDisplay::resetProjection()
{
	observationSettings().focus[0] =  0.0f;
	observationSettings().focus[1] = 0.0f;
	observationSettings().focus[2] = 0.0f;
	observationSettings().rotation[0] = 0.0f;
	observationSettings().rotation[1] = 0.0f;
	observationSettings().rotation[2] = 0.0f;
	observationSettings().rotation[3] = 1.0f;
	if ((observationSettings().projection[0]!=observationSettings().projection[1]) &&
		(observationSettings().projection[1]!=observationSettings().projection[2]) &&
		(observationSettings().projection[2]!=observationSettings().projection[0]))
	{
		idx_type r = this->input.points.rows();
		this->projection = Matrix::Zero(3, r);
		for (idx_type i=0; i<3&&i<r; ++i)
		{
			idx_type index = observationSettings().projection[i];
			if (index>0 && index-1<r)
				this->projection(i,index-1)=1;
		}
	}
	adjustDisplay();
}

Matrix GMMDisplay::randomOrthogonalProjection(idx_type from, idx_type to)
{
	if (from<=to)
		return Matrix::Identity(to, from);
		
	std::mt19937 gen(time(0));

	std::normal_distribution<> dis(0, 1);
	Matrix m(from, to);

	for (idx_type i=0; i<from; ++i)
		for (idx_type j=0; j<to; ++j)
			m(i,j)=dis(gen);

	linalg::gramschmidt(m);
	return m.transpose();
}

void GMMDisplay::draw()
{  
	double mat[4*4]; // rotation matrix
	idx_type n = this->input.points.cols();
	idx_type d = this->input.points.rows();

	if(observationSettings().alpha != this->alpha)
	{
		this->alpha = observationSettings().alpha;
		if(this->alpha > 0)
			this->quantile = boost::math::quantile(this->csd,this->alpha);
		else
			this->quantile = 1.;
	}
		
	if (observationSettings().random)
	{
		observationSettings().random = false;
		this->projection = randomOrthogonalProjection(d,3);
		adjustDisplay();
	}

	if (observationSettings().reset)
	{
		observationSettings().reset = false;
		resetProjection();
	}

	if (observationSettings().flip)
	{
		observationSettings().flip = false;
		for (size_t i=0; i<d; ++i)
			projection(0,i)*=-1;
		adjustDisplay();
	}
	
	if(this->pointSize != displaySettings().pointSize)
	{
		adjustDisplay();
		this->pointSize = displaySettings().pointSize;
	}

	if (observationSettings().write)
	{
		observationSettings().write = false;
		int g = observationSettings().component;
		g--;
			GMMIOUtil::appendToGMM(gmmlabSettings().outputFile,
				this->solutions[observationSettings().solution],
				this->tags[observationSettings().solution]);
			if (commonSettings().verbose)
			{
				std::ostringstream stream;
				stream << "solution " << observationSettings().solution;
				if (g>=0)
					stream << "." << g;
			}
	}

	// Clear frame buffer using bgColor
	glClearColor(displaySettings().bgColor[0], displaySettings().bgColor[1], displaySettings().bgColor[2], 1);
	glClear( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(displaySettings().fovy, displaySettings().aspectRatio, displaySettings().zNear, displaySettings().zFar);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if (observationSettings().follow
		&& observationSettings().component > 0
		&& observationSettings().solution < this->solutions.size()
		&& observationSettings().component-1 < this->solutions[observationSettings().solution].means.cols())
	{
		Vector p = this->projection*this->solutions[observationSettings().solution].means.col(
			observationSettings().component-1);
		observationSettings().focus[0] = p[0];
		observationSettings().focus[1] = p[1];
		observationSettings().focus[2] = p[2];
	}

	gluLookAt(0,0,pow(2,observationSettings().zoomLog)*displaySettings().distance,
    			0,0,0,
    			0,1,0);

	convertQuaternionToMatrix(observationSettings().rotation, mat);
	glMultMatrixd(mat);

	glTranslated(-observationSettings().focus[0], -observationSettings().focus[1],
		-observationSettings().focus[2]);

	// draw	
	if (!observationSettings().hideSamples && observationSettings().solution < this->sampleSets.size())
		drawBoxes(this->sampleSets[observationSettings().solution]);

	if (observationSettings().solution < this->solutions.size())
	{
		if (!observationSettings().hideSolution)
			drawGMM(this->solutions[observationSettings().solution],
				displaySettings().ellipsoidColor, displaySettings().highlightColor);		
		else if (observationSettings().component > 0
				&& observationSettings().component-1 < this->solutions[observationSettings().solution].weights.size())
			drawGaussian(this->solutions[observationSettings().solution].means.col(observationSettings().component-1),
				this->solutions[observationSettings().solution].covariances[observationSettings().component-1],
				displaySettings().highlightColor);
	}

	if (observationSettings().showTruth && !this->truth.empty())
		drawGMM(this->truth, displaySettings().truthColor);

	if (!observationSettings().hideInput)
	{
		
		if (this->pointsDL!=0 && !observationSettings().partition)
		{
			glColor4ubv(displaySettings().pointColor);
			glCallList(this->pointsDL);
		}
		else
			drawPoints(this->input.points);
	}
	
}

void GMMDisplay::initGL()
{
	
     // ----
     // (1) solid
     glEnable(GL_DEPTH_TEST); // enables culling
     
     // (2) transparent ellipses
//   glDisable(GL_CULL_FACE);
//   glEnable(GL_BLEND);
//   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
     // ----
    
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2.0);
    glPointSize(displaySettings().pointSize);

    // Set light
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);    // use default light diffuse and position
    glEnable(GL_NORMALIZE);
    float v[4]; // will be used to set light parameters
    v[0] = v[1] = v[2] = .7f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, v);
    v[0] = v[1] = v[2] = .1f; v[3] = 1.0f;
    glLightfv(GL_LIGHT0, GL_AMBIENT, v);
}



int GMMDisplay::initGLFW()
{
    const GLFWvidmode* mode;   // GLFW video mode

    // Intialize GLFW
    if( !glfwInit() )
    {
        // An error occured
        fprintf(stderr, "GLFW initialization failed\n");
        return 1;
    }

    // Create a window
    mode = glfwGetVideoMode(glfwGetPrimaryMonitor());    
    int f = displaySettings().fullscreen?1:2;
    int width = (int)(mode->width/f);
    int height =  (int)(mode->height/f);
    window = glfwCreateWindow(width, height, "GMMLab", NULL, NULL); //displaySettings().fullscreen?glfwGetPrimaryMonitor():NULL, NULL);
    if( !window )
	   //glfwOpenWindow(mode->width/f, mode->height/f, mode->redBits, mode->greenBits, mode.BlueBits,
           //             0, 16, 0, displaySettings().fullscreen?GLFW_FULLSCREEN:GLFW_WINDOW ) )
    {
        // A fatal error occured
        fprintf(stderr, "Cannot open GLFW window\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwWindowHint(GLFW_RED_BITS, mode->redBits);
    glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
    glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
    glfwWindowHint(GLFW_ALPHA_BITS, 0);
    glfwWindowHint(GLFW_DEPTH_BITS, 16);
    glfwWindowHint(GLFW_STENCIL_BITS, 0);
    
    
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_FALSE); //   glfwEnable(GLFW_KEY_REPEAT);

    // Initialize AntTweakBar
    TwInit(TW_OPENGL, NULL);
    WindowSizeCB(window, width, height);

    // Add Message to the help bar.
	TwDefine(" GLOBAL help='This Tool shows the solutions computed by the Algorithm.' ");
	TwDefine(" GLOBAL iconpos=tl ");
//	TwDefine(" GLOBAL fontsize=3 "); // use large font - SEEMS BROKEN
//	TwDefine(" observation size='200 400' color='96 216 224' ");

	displayBar = new DisplayBar(displaySettings());
	observationBar = new ObservationBar(observationSettings());

	std::stringstream opt_cost_label;
	opt_cost_label << " label='opt-cost (" << commonSettings().costmeasure << ")' help='Cost of underlying truth (if available).' ";
	std::stringstream cost_label;
	cost_label << "label='cost (" << commonSettings().costmeasure << ")' help='Cost of selected solution.' ";
	
	TwBar* bar = TwNewBar("data");
	TwDefine(" data label='Statistics' position='40 40' size='490 140' valueswidth=300 iconified=false ");
	TwAddVarRO(bar, "tag", TW_TYPE_STDSTRING, &this->statistics.tag,
				" label='tag' help='Tag of selected solution.' ");
	TwAddVarRO(bar, "time", TW_TYPE_DOUBLE, &this->statistics.runtime,
				" label='comp.time' help='Computation time of selected solution.' ");
	TwAddVarRO(bar, "nll", TW_TYPE_DOUBLE, &this->statistics.nll,
				cost_label.str().c_str());
	TwAddSeparator(bar, NULL, NULL);
	TwAddVarRO(bar, "truthnnl", TW_TYPE_DOUBLE, &this->statistics.truthNLL,
				opt_cost_label.str().c_str());
	TwAddVarRO(bar, "truthsep", TW_TYPE_DOUBLE, &this->statistics.separation,
				" label='gmm-separation' help='Separation of underlying GMM (if available).' ");

    // Set GLFW event callbacks
    // - Redirect window size changes to the callback function WindowSizeCB
    glfwSetWindowSizeCallback(window, WindowSizeCB);
    
   
    // - Directly redirect GLFW mouse button events to AntTweakBar
    glfwSetMouseButtonCallback(window, (GLFWmousebuttonfun)(twEventMouseButtonGLFW3));
    // - Directly redirect GLFW mouse position events to AntTweakBar
    glfwSetCursorPosCallback(window, (GLFWcursorposfun)(twEventMousePosGLFW3));
    // - Directly redirect GLFW mouse wheel events to AntTweakBar
//    glfwSetMouseWheelCallback((GLFWmousewheelfun)TwEventMouseWheelGLFW);
    // - Directly redirect GLFW key events to AntTweakBar
    glfwSetKeyCallback(window, (GLFWkeyfun)(twEventKeyGLFW3));
    // - Directly redirect GLFW char events to AntTweakBar
    glfwSetCharCallback(window, (GLFWcharfun)(twEventCharGLFW3));
}


double calcFPS(GLFWwindow* window, double theTimeInterval = 1.0, std::string theWindowTitle = "NONE")
{
	// Static values which only get initialised the first time the function runs
	static double t0Value       = glfwGetTime(); // Set the initial time to now
	static int    fpsFrameCount = 0;             // Set the initial FPS frame count to 0
	static double fps           = 0.0;           // Set the initial FPS value to 0.0
 
	// Get the current time in seconds since the program started (non-static, so executed every time)
	double currentTime = glfwGetTime();
 
	// Ensure the time interval between FPS checks is sane (low cap = 0.1s, high-cap = 10.0s)
	// Negative numbers are invalid, 10 fps checks per second at most, 1 every 10 secs at least.
	if (theTimeInterval < 0.1)
	{
		theTimeInterval = 0.1;
	}
	if (theTimeInterval > 10.0)
	{
		theTimeInterval = 10.0;
	}
 
	// Calculate and display the FPS every specified time interval
	if ((currentTime - t0Value) > theTimeInterval)
	{
		// Calculate the FPS as the number of frames divided by the interval in seconds
		fps = (double)fpsFrameCount / (currentTime - t0Value);
 
		// If the user specified a window title to append the FPS value to...
		if (theWindowTitle != "NONE")
		{
			// Convert the fps value into a string using an output stringstream
			std::ostringstream stream;
			stream << fps;
			std::string fpsString = stream.str();
 
			// Append the FPS value to the window title details
			theWindowTitle += " | FPS: " + fpsString;
 
			// Convert the new window title to a c_str and set it
			const char* pszConstString = theWindowTitle.c_str();
			glfwSetWindowTitle(window, pszConstString);
		}
		else // If the user didn't specify a window to append the FPS to then output the FPS to the console
		{
			std::cout << "FPS: " << fps << std::endl;
		}
 
		// Reset the FPS frame counter and set the initial time to be now
		fpsFrameCount = 0;
		t0Value = glfwGetTime();
	}
	else // FPS calculation time interval hasn't elapsed yet? Simply increment the FPS frame counter
	{
		fpsFrameCount++;
	}
 
	// Return the current FPS - doesn't have to be used if you don't want it!
	return fps;
}

void GMMDisplay::mainloop()
{  
	// Main loop (repeated while window is not closed and [ESC] is not pressed)
	while( !glfwWindowShouldClose(window) ) // !glfwGetKey(window, GLFW_KEY_ESCAPE) &&
	{  
		// Draw solution
		draw();

		if (observationSettings().solution < solutions.size())
		{
			observationBar->setComponentRange(0,this->solutions[observationSettings().solution].weights.size());

			this->statistics.tag = observationSettings().solution < this->tags.size()?
				this->tags[observationSettings().solution]:"";
			this->statistics.nll = observationSettings().solution < this->nlls.size()?
				this->nlls[observationSettings().solution]:-1;
			this->statistics.runtime = observationSettings().solution < this->runtimes.size()?
				this->runtimes[observationSettings().solution]:-1;
		}

		// Draw tweak bars
		TwDraw();

		// Present frame buffer
		glfwSwapBuffers(window);
		glfwPollEvents();
		
		// update FPS
		calcFPS(window, 1, "GMM Lab");
	}
}