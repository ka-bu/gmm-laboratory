#include "../base.h"

#ifndef GMMDISPLAY_H
#define GMMDISPLAY_H

#include <boost/math/distributions/chi_squared.hpp>

#include "../base/parameters.h"
#include "../base/commonutil.h"

#include "../settings/settings.h"
#include "../settings/displaysettings.h"
#include "../settings/observationsettings.h"

#include "displaybar.h"
#include "observationbar.h"

#ifndef GLFW_INCLUDE_GLU
#define GLFW_INCLUDE_GLU
#endif
#include <GLFW/glfw3.h>

/**
* @brief Display
*/
class GMMDisplay
{
public:
	GMMDisplay(Parameters& truth, commonutil::DataSet const&, std::vector<commonutil::DataSet> const&,
		std::vector<Parameters> const&, std::vector<std::string> const&, std::vector<double> const&,
		std::vector<fp_type> const&);
	
	virtual ~GMMDisplay();

	void mainloop();

private:
	static void WindowSizeCB(GLFWwindow*, int, int);
	void prepareDisplayLists();
	GLdouble* convertMatrix(Matrix const&);
	void drawGaussian(Vector const& mean, Matrix covariance, float const* color);
	void drawGMM(Parameters const&, float const*, float const* = NULL);
	void drawPoints(Matrix const& points);
	void drawBoxes(commonutil::DataSet const&);
	void draw();
	void convertQuaternionToMatrix(const double*, double*);
	void adjustDisplay();
	void resetProjection();
	Matrix randomOrthogonalProjection(idx_type, idx_type);
	void initGL();
	int initGLFW();
	
	Parameters& truth;
	commonutil::DataSet const& input;
	std::vector<commonutil::DataSet> const& sampleSets;
	std::vector<Parameters> const& solutions;
	std::vector<std::string> const& tags;
	std::vector<fp_type> const& nlls;
	std::vector<double> const& runtimes;
	fp_type quantile;
	fp_type alpha;
	
	boost::math::chi_squared_distribution<fp_type> csd;
		
	GLFWwindow* window;
	DisplayBar* displayBar;
	ObservationBar* observationBar;

	Matrix projection;
	GLuint pointsDL, hiSphereDL, loSphereDL, boxDL;
	fp_type lightDistance;
	double pointSize;
	
	struct {
		fp_type truthNLL;
		std::string tag;
		fp_type nll;
		double runtime;
		fp_type separation;
	} statistics;
	
	inline static void twEventMouseButtonGLFW3(GLFWwindow* window, int button, int action, int mods)
	{
		TwEventMouseButtonGLFW(button, action);
	}

	inline static void twEventMousePosGLFW3(GLFWwindow* window, double xpos, double ypos)
	{
		TwEventMousePosGLFW((int)xpos, (int)ypos);
	}

	inline static void twEventKeyGLFW3(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		TwEventKeyGLFW(key, action);
	}

	inline static void twEventCharGLFW3(GLFWwindow* window, unsigned int codepoint)
	{
		TwEventCharGLFW(codepoint, GLFW_PRESS);
}
	
};

#endif
