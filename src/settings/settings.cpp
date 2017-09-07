#include "settings.h"

CommonSettings& commonSettings()
{
	static CommonSettings cs;
	return cs;
}

GMMLabSettings& gmmlabSettings()
{
	static GMMLabSettings gs;
	return gs;
}

TestLabSettings& testlabSettings()
{
	static TestLabSettings ts;
	return ts;
}

DisplaySettings& displaySettings()
{
	static DisplaySettings ds;
	return ds;
}

ObservationSettings& observationSettings()
{
	static ObservationSettings os;
	return os;
}

