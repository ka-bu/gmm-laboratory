#include "../base.h"

#include "displaybar.h"

DisplayBar::DisplayBar(DisplaySettings& settings)
{
	this->bar = TwNewBar("display");
	
	TwAddVarRW(this->bar, "bgColor", TW_TYPE_COLOR3F, &settings.bgColor, " label='BG color' ");

	TwAddVarRW(this->bar, "pointColor", TW_TYPE_COLOR32, &settings.pointColor, 
       			" label='Point color' alpha help='Color and transparency of the points.' ");

	TwAddVarRW(this->bar, "boxColor", TW_TYPE_COLOR3F, &settings.boxColor, 
       			" label='Box color' help='Color and transparency of the boxes.' ");

	TwAddVarRW(this->bar, "truthColor", TW_TYPE_COLOR4F, &settings.truthColor, 
				" label='Truth color' help='Color and transparency of the underlying truth.' ");

	TwAddVarRW(this->bar, "ellipsoidColor", TW_TYPE_COLOR3F, &settings.ellipsoidColor, 
				" label='Gaussian color' help='Color and transparency of the ellipsoids.' ");

	TwAddVarRW(this->bar, "selectedColor", TW_TYPE_COLOR3F, &settings.highlightColor, 
				" label='Selection color' help='Color and transparency of the selected ellipsoid.' ");

	TwAddVarRW(this->bar, "wire", TW_TYPE_BOOLCPP, &settings.wireframe, 
				" label='Wireframe' key=W help='Toggle wireframe display mode.' ");

	TwAddVarRW(this->bar, "boxlength", TW_TYPE_DOUBLE, &settings.length, 
				" label='Box size' min=-10 max=10 step=0.1 keyIncr=b keyDecr=B ");
	
	TwAddVarRW(this->bar, "pointSize", TW_TYPE_DOUBLE, &settings.pointSize, " label='Point size' min=1 max=1000 step=1 keyIncr=. keyDecr=: ");

	TwAddVarRW(this->bar, "fovy", TW_TYPE_DOUBLE, &settings.fovy, 
				" label='view angle' group=perspective ");

	TwAddVarRW(this->bar, "aspect", TW_TYPE_DOUBLE, &settings.aspectRatio, 
				" label='aspect ratio' group=perspective ");

	TwAddVarRW(this->bar, "znear", TW_TYPE_DOUBLE, &settings.zNear, 
				" label='z-Near' group=perspective keyIncr=N keyDecr=n ");

	TwAddVarRW(this->bar, "zfar", TW_TYPE_DOUBLE, &settings.zFar, 
				" label='z-Far' group=perspective keyIncr=m keyDecr=M ");

	TwDefine(" display/perspective label='Perspective' "); 
         
	TwDefine(" display label='Display settings' position='40 600' size='200 280' iconified=false ");
}