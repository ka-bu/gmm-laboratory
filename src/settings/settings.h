#ifndef SETTINGS_H
#define SETTINGS_H

#include "commonsettings.h"
#include "gmmlabsettings.h"
#include "testlabsettings.h"
#include "displaysettings.h"
#include "observationsettings.h"

CommonSettings& commonSettings();
GMMLabSettings& gmmlabSettings();
TestLabSettings& testlabSettings();
DisplaySettings& displaySettings();
ObservationSettings& observationSettings();

#endif