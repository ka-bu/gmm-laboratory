#ifndef OBSERVATIONSETTINGS_H
#define OBSERVATIONSETTINGS_H

struct ObservationSettings
{
	unsigned int solution; // number of displayed solution
	unsigned int component; // number of highlighted component
	bool hideSolution; // hide the ellipsoids of the solutions
	bool hideSamples; // hide the samples
	bool hideInput; // hide the input
	bool showTruth; // show the underlying truth if not NULL ?
	bool follow; //  should the focus follow the selection ?
	double focus[3]; // coordinates of the point in focus
	double rotation[4]; // orientation (stored as quaternion)
	double zoomLog; // logarithm of zoom factor
	unsigned int projection[3]; // indices for resetting the projection
	bool reset; // 'Reset Projection' button pressed
	bool random; // 'Random Projection' button pressed
	bool flip; // inside out
	bool write; // inside out
	double alpha;
	double meansOnly; // set to 0 if *not* only means shall be displayed; size of the sphere around the mean otherwise
	bool partition;

	ObservationSettings()
	{
		solution = 0;	// number of displayed solution
		component = 0;	// number of highlighted gaussian
		hideSolution = false;	// hide the ellipsoids of the solutions
		hideSamples = false;	// hide the samples
		hideInput = false;		// hide the input
		showTruth = true;		// show the underlying truth if not NULL ?
		follow = false;		// should the focus follow the selection ?
		focus[0] =  0.0f;
		focus[1] = 0.0f;
		focus[2] = 0.0f;			// coordinates of the point in focus
		rotation[0] = 0.0f;
		rotation[1] = 0.0f;
		rotation[2] = 0.0f;
		rotation[3] = 1.0f; 	// orientation (stored as quaternion)
		zoomLog = 0;			// logarithm of zoom factor
		projection[0] = 1;
		projection[1] = 2;
		projection[2] = 3;
		reset = false;			// 'Reset Projection' button pressed
		random = false;		// 'Random Projection' button pressed
		flip = false;			// inside out
		write = false;			// inside out
		alpha = 0.0f;
		meansOnly = 0.0f;
		partition = false;
	}
};

#endif
