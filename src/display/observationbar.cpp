#include "../base.h"

#include "observationbar.h"

ObservationBar::ObservationBar(ObservationSettings& settings)
{
	this->bar = TwNewBar("observation");

	TwAddVarRW(this->bar, "solution", TW_TYPE_INT32, &settings.solution,
				" label='Solution' keyIncr=s keyDecr=S help='Changes the displayed solution.' ");

	TwAddVarRW(this->bar, "component", TW_TYPE_INT32, &settings.component,
				" label='Component' keyIncr=g keyDecr=G help='Changes the selected component.' ");

	TwAddVarRW(this->bar, "follow", TW_TYPE_BOOLCPP, &settings.follow,
				" label='Follow selection' key=f help='Let the focus follow the selection.' ");

	TwAddVarRW(this->bar, "hide_solution", TW_TYPE_BOOLCPP, &settings.hideSolution,
				" label='Hide solution' key=h help='Hide the ellipsoids of the solution.' ");

	TwAddVarRW(this->bar, "hide_samples", TW_TYPE_BOOLCPP, &settings.hideSamples,
				" label='Hide samples' key=H help='Hide the samples drawn by the algorithm.' ");

	TwAddVarRW(this->bar, "hide_input", TW_TYPE_BOOLCPP, &settings.hideInput,
				" label='Hide input' key=i help='Hide the input points.' ");

	TwAddVarRW(this->bar, "truth", TW_TYPE_BOOLCPP, &settings.showTruth,
				" label='Show truth' key=t help='Show the underlying truth if available.' ");

	TwAddVarRW(this->bar, "zoom", TW_TYPE_DOUBLE, &settings.zoomLog,
				" label='Distance' min=-10 max=10 step=0.1 keyIncr=D keyDecr=d ");

	TwAddVarRW(this->bar, "alpha", TW_TYPE_DOUBLE, &settings.alpha, " label='alpha' min=0 max=1 step=0.01 help='show (1-alpha)*100%-ellipsoids'");
		
	TwAddVarRW(this->bar, "means_only", TW_TYPE_DOUBLE, &settings.meansOnly, " label='means_only' min=0 max=10000 step=50 help='show spheres of constant size around means'");
	
	TwAddVarRW(this->bar, "partition", TW_TYPE_BOOLCPP, &settings.partition, " label='partition' help='show spheres of constant size around means'");
	
	TwAddVarRW(this->bar, "rotation", TW_TYPE_QUAT4D, &settings.rotation,
				" label='Rotation' opened=true help='Change the object orientation.' ");

	TwAddVarRW(this->bar, "focusx", TW_TYPE_DOUBLE, &settings.focus[0],
				" label='X' group=focus step=0.1 keyIncr=x keyDecr=X help='Where to look at.' ");

	TwAddVarRW(this->bar, "focusy", TW_TYPE_DOUBLE, &settings.focus[1],
				" label='Y' group=focus step=0.1 keyIncr=y keyDecr=Y help='Where to look at.' ");

	TwAddVarRW(this->bar, "focusz", TW_TYPE_DOUBLE, &settings.focus[2],
				" label='Z' group=focus step=0.1 keyIncr=z keyDecr=Z help='Where to look at.' ");

	TwDefine(" observation/focus label='Looking at' opened=false ");

	TwAddVarRW(this->bar, "first", TW_TYPE_INT32, &settings.projection[0],
				" label='X' group=projection step=1 keyIncr=2 keyDecr=1 help='Coordinate for projection onto X' ");

	TwAddVarRW(this->bar, "second", TW_TYPE_INT32, &settings.projection[1],
				" label='Y' group=projection step=1 keyIncr=4 keyDecr=3 help='Coordinate for projection onto Y' ");

	TwAddVarRW(this->bar, "third", TW_TYPE_INT32, &settings.projection[2],
				" label='Z' group=projection step=1 keyIncr=6 keyDecr=5 help='Coordinate for projection onto Z' ");

	TwAddButton(this->bar, "reset", ObservationBar::buttonCallback, &settings.reset,
                " label='Reset' group=projection key=r help='reset projection.' ");

	TwAddButton(this->bar, "random", ObservationBar::buttonCallback, &settings.random,
                " label='Create randomly' group=projection key=p help='generate random projection.' ");

	TwAddButton(this->bar, "flip", ObservationBar::buttonCallback, &settings.flip,
                " label='Flip inside out' group=projection key=q help='flip orientation of the projection.' ");

	TwAddButton(this->bar, "write", ObservationBar::buttonCallback, &settings.write,
                " label='Write solution' key=w help='Write solution to file.' ");

	TwDefine(" observation/projection label='3D projection' ");

	TwDefine(" observation label='Observation settings' position='40 200' size='250 370' valueswidth=100 iconified=false ");
}

void TW_CALL ObservationBar::buttonCallback(void *clientData)
{
    bool* b = static_cast<bool*>(clientData);
    *b = true;
}

void ObservationBar::setSolutionRange(int min, int max)
{
	TwSetParam(this->bar, "solution", "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(this->bar, "solution", "max", TW_PARAM_INT32, 1, &max);
}

void ObservationBar::setComponentRange(int min, int max)
{
	TwSetParam(this->bar, "component", "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(this->bar, "component", "max", TW_PARAM_INT32, 1, &max);
}

void ObservationBar::setProjectionRange(int min, int max)
{
	TwSetParam(this->bar, "first", "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(this->bar, "first", "max", TW_PARAM_INT32, 1, &max);

	TwSetParam(this->bar, "second", "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(this->bar, "second", "max", TW_PARAM_INT32, 1, &max);

	TwSetParam(this->bar, "third", "min", TW_PARAM_INT32, 1, &min);
	TwSetParam(this->bar, "third", "max", TW_PARAM_INT32, 1, &max);
}