#ifndef DISPLAYSETTINGS_H
#define DISPLAYSETTINGS_H

struct DisplaySettings
{
	float bgColor[3];				// background color 
	unsigned char pointColor[4]; 	// point color (32bits RGBA)
	float boxColor[4]; 				// mean color (32bits RGBA)
	float ellipsoidColor[4];		// ellipsoid color (32bits RGBA)
	float highlightColor[4];		// highlight color (32bits RGBA)
	float truthColor[4];			// truth color (32bits RGBA)
	bool wireframe;					// only draw wireframe?
	double length; 
	double distance; 
	double fovy; 
	double aspectRatio;
	double zNear; 
	double zFar;
	double pointSize;
	bool fullscreen; 
	bool coloredCluster = false;
	float clusterColors[24][4];

	DisplaySettings()
	{
/*

// option 1: 

 		bgColor[0] = 0;	// background color 
 		bgColor[1] = 0;
 		bgColor[2] = .3;


		pointColor[0] = 255; // point color (32bits RGBA)
		pointColor[1] = 255;
		pointColor[2] = 255;
		pointColor[3] = 64;

		boxColor[0] = 0; // mean color (32bits RGBA)
		boxColor[1] = 1;
		boxColor[2] = 1;
		boxColor[3] = 1;
		
		ellipsoidColor[0] = 1; // ellipsoid color (32bits RGBA)
		ellipsoidColor[1] = 0;
		ellipsoidColor[2] = 0;
		ellipsoidColor[3] = 1;
		
// end of option 1
*/

// option 2: 

		bgColor[0] = 255;	// background color 
 		bgColor[1] = 255;
 		bgColor[2] = 255;		
		
		pointSize = 2;
		
		pointColor[0] = 0; // point color (32bits RGBA)
		pointColor[1] = 0;
		pointColor[2] = 0;
		pointColor[3] = 120;
		
		boxColor[0] = 0; // mean color (32bits RGBA)
		boxColor[1] = 0;
		boxColor[2] = 0;
		boxColor[3] = 255;		
		
		ellipsoidColor[0] = 1; // ellipsoid color (32bits RGBA)
		ellipsoidColor[1] = 0;
		ellipsoidColor[2] = 0;
		ellipsoidColor[3] = 50;
		
// end of option 2
		
		highlightColor[0] = 1; // highlight color (32bits RGBA)
		highlightColor[1] = .7;
		highlightColor[2] = 0;
		highlightColor[3] = 1;
		
		truthColor[0] = 0; // truth color (32bits RGBA)
		truthColor[1] = 1;
		truthColor[2] = 0;
		truthColor[3] = 0.2; 
		
		fp_type colvalue = 1;
		for (idx_type i = 0; i < 4; ++i)
		{
			for (idx_type j = 0; j < 6; ++j)
			{
				clusterColors[i*6+j][0] = 0;
				clusterColors[i*6+j][1] = 0;
				clusterColors[i*6+j][2] = 0;
				switch(j)
				{
					case 0:
						clusterColors[i*6+j][0] = colvalue;
						break;
					case 1:
						clusterColors[i*6+j][1] = colvalue;
						break;
					case 2:
						clusterColors[i*6+j][2] = colvalue;
						break;
					case 3:
						clusterColors[i*6+j][0] = colvalue;
						clusterColors[i*6+j][1] = colvalue;
						break;
					case 4:
						clusterColors[i*6+j][1] = colvalue;
						clusterColors[i*6+j][2] = colvalue;
						break;
					case 5:
						clusterColors[i*6+j][0] = colvalue;
						clusterColors[i*6+j][2] = colvalue;
						break;
				}
				clusterColors[i*6+j][3] = 1;
			}
			colvalue *= 0.5;
		}
		
		wireframe = false;	// only draw wireframe?
		length = -3;
		distance = 1;
		fovy = 40;
		aspectRatio = 1;
		zNear = 1;
		zFar = 2000;
		fullscreen = false;
	}
};

#endif
