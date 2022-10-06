#include "SBG.hpp"

using namespace time_literals;

#define EKF2_MAX_INSTANCES 9

SBG::SBG() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::nav_and_controllers)
{}

// Static callback
SbgErrorCode SBG::onLogReceived(SbgEComHandle *pHandle, SbgEComClass msgClass, SbgEComMsgId msg, const SbgBinaryLogData *pLogData, void *pUserArg)
{
	if (!SBG::instance) {
		return SBG_NULL_POINTER;
	}
	SBG::instance->onLog(msg, pLogData);
	return SBG_NO_ERROR;
}

bool SBG::init()
{
	// Connect to module
	int ret = connect();
	if (ret) {
		PX4_ERR("[SBG] Connection error");
		return;
	}

	// Raw sensors
	_sbg_ekf_nav_pub.advertise();
	_sbg_ekf_quat_pub.advertise();

	// Estimator
	_estimator_selector_status_pub.advertise();
	_sensor_selection_pub.advertise();
	_vehicle_attitude_pub.advertise();
	_vehicle_global_position_pub.advertise();
	_vehicle_local_position_pub.advertise();
	_vehicle_odometry_pub.advertise();
	_wind_pub.advertise();

	// if (!_trigger_sub.registerCallback()) {
	// 	PX4_ERR("callback registration failed");
	// 	return false;
	// }


	SbgErrorCode errorCode = SBG_NO_ERROR;
	// printf("[SBG INTERFACE] %p\n", &_interface);
	// errorCode = sbgInterfaceSerialCreate(&_interface, "/dev/tty", 9600);
	printf("[SBG CREATED] %p %p %p; error: %d\n", &_comHandle, &_interface, &_deviceInfo, errorCode);
	ScheduleOnInterval(500_ms);


	return true;
}

int SBG::connect()
{
	SbgErrorCode errorCode;
	char* port = "/dev/tty0";
	errorCode = sbgInterfaceSerialCreate(&_interface, port, 9600);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to create the interface");
		return PX4_ERROR;
	}

	// Create the sbgECom library and associate it with the created interfaces
	errorCode = sbgEComInit(&_comHandle, &_interface);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to initialize the sbgECom library");
		return PX4_ERROR;
	}

	// Define callbacks for received data
	sbgEComSetReceiveLogCallback(&_comHandle, SBG::onLogReceived, NULL);

	// Configure output logs to 200 Hz
	errorCode = sbgEComCmdOutputSetConf(&_comHandle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_EKF_QUAT, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_EKF_QUAT");
		return PX4_ERROR;
	}

	errorCode = sbgEComCmdOutputSetConf(&_comHandle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_IMU_DATA, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_IMU_DATA");
		return PX4_ERROR;
	}

	errorCode = sbgEComCmdOutputSetConf(&_comHandle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_EKF_NAV, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_EKF_NAV");
		return PX4_ERROR;
	}

	// Display device information
	errorCode = sbgEComCmdGetInfo(&_comHandle, &_deviceInfo);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to get device information");
		return PX4_ERROR;
	}
	PX4_INFO("[SBG] Connected. SN: %.9u, version: %s\n", _deviceInfo.serialNumber, SBG_E_COM_VERSION_STR);

	return PX4_OK;
}


void SBG::onLog(const SbgEComMsgId &msg, const SbgBinaryLogData *pLogData)
{
	switch (msg)
	{
	case SBG_ECOM_LOG_EKF_EULER:
	{
		for (int i = 0; i < 3; i++)
		{
			pLogData->ekfEulerData.eulerStdDev[i];
		}
		break;
	}
	case SBG_ECOM_LOG_EKF_QUAT:
	{
		pLogData->ekfQuatData.quaternion[1];
		pLogData->ekfQuatData.quaternion[2];
		pLogData->ekfQuatData.quaternion[3];
		pLogData->ekfQuatData.quaternion[0];
		break;
	}
	case SBG_ECOM_LOG_IMU_DATA:
	{
		for (int i = 0; i < 3; i++)
		{
			pLogData->imuData.gyroscopes[i];
			pLogData->imuData.accelerometers[i];
		}
		break;
	}

	case SBG_ECOM_LOG_GPS1_HDT:
	{
		// _data.gnssHdtValid = false;
		if (sbgEComLogGpsHdtGetStatus(pLogData->gpsHdtData.status) == SBG_ECOM_HDT_SOL_COMPUTED)
		{
			// _data.gnssHdt = pLogData->gpsHdtData.heading * D2R - M_PI; // M_PI is yaw offset between baseline and beam.
			// _data.gnssHdtAcc = pLogData->gpsHdtData.headingAccuracy * D2R;
			// _data.gnssHdtValid = true;
		}
		break;
	}
	case SBG_ECOM_LOG_GPS1_POS:
	{
		const SbgLogGpsPos &gpsData = pLogData->gpsPosData;
		if(sbgEComLogGpsPosGetStatus(gpsData.status) == SBG_ECOM_POS_SOL_COMPUTED)
		{
			// _data.gnssTimestamp = gpsData.timeOfWeek; // ~5Hz only
			// _data.gnssCoord = gnss::Coord(gpsData.latitude, gpsData.longitude,	gpsData.altitude);
			// _data.gnssCoordValid = true;
		}
		break;
	}
	case SBG_ECOM_LOG_EKF_NAV:
	{
		pLogData->ekfNavData.position;
		pLogData->ekfNavData.positionStdDev;
		pLogData->ekfNavData.velocity;
		pLogData->ekfNavData.velocityStdDev;
		break;
	}
	default:
		break;
	}
}

void SBG::publishEstimatorSelectorStatus()
{
	estimator_selector_status_s selector_status{};
	selector_status.primary_instance = 0;
	selector_status.instances_available = 1;
	selector_status.instance_changed_count = 0;
	selector_status.last_instance_change = 0;
	selector_status.accel_device_id = 0;
	selector_status.baro_device_id = 0;
	selector_status.gyro_device_id = 0;
	selector_status.mag_device_id = 0;
	selector_status.gyro_fault_detected = 0;
	selector_status.accel_fault_detected = 0;

	for (int i = 0; i < DIM(selector_status.healthy); i++) {
		selector_status.combined_test_ratio[i] = 0;
		selector_status.relative_test_ratio[i] = 0;
		selector_status.healthy[i] = true;
	}

	for (int i = 0; i < DIM(selector_status.accumulated_gyro_error); i++) {
		selector_status.accumulated_gyro_error[i] = 0;
		selector_status.accumulated_accel_error[i] = 0;
	}

	selector_status.timestamp = hrt_absolute_time();
	_estimator_selector_status_pub.publish(selector_status);
}

void SBG::publishVehicleAttitude()
{
	vehicle_attitude_s attitude{};
	attitude.q[0] = 0;
	attitude.q[1] = 0;
	attitude.q[2] = 0;
	attitude.q[3] = 0;
	//
	attitude.delta_q_reset[0] = 0;
	attitude.delta_q_reset[1] = 0;
	attitude.delta_q_reset[2] = 0;
	attitude.delta_q_reset[3] = 0;
	attitude.quat_reset_counter = 0;

	attitude.timestamp = hrt_absolute_time();
}

void SBG::publishVehicleLocalPosition()
{
	vehicle_local_position_s local_position{};

	local_position.timestamp = hrt_absolute_time();
}
void SBG::publishVehicleGlobalPosition()
{
	vehicle_global_position_s global_position{};
	// global_position.alt_ellipsoid

	global_position.timestamp = hrt_absolute_time();
}

void SBG::publishVehicleOdometry()
{
	vehicle_odometry_s odometry{};
	odometry.timestamp = hrt_absolute_time();
}


void SBG::Run()
{
	if (should_exit()) {
		// _trigger_sub.unregisterCallback();
		exit_and_cleanup();
		return;
	}

	//
	// px4_usleep(500000);
	printf("[SBG loop]\n");
	sbg_ekf_nav_s nav{};
	nav.position[0] = 1.;
	nav.position[1] = 2.;
	nav.position[2] = 3.;
	_sbg_ekf_nav_pub.publish(nav);

	// Publish estimator status
	estimator_selector_status_s selector_status{};


}


int SBG::task_spawn(int argc, char *argv[])
{
	if (SBG::instance) {
		PX4_ERR("SBG task already spawned");
		return PX4_ERROR;
	}
	SBG::instance = new SBG();
	if (!SBG::instance) {
		PX4_ERR("SBG allocation failed");
		return PX4_ERROR;
	}
	SBG::instance->init();
	_object.store(SBG::instance);
	_task_id = task_id_is_work_queue;
	return PX4_OK;
}

int SBG::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int SBG::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description


)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("sbg", "system");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

extern "C" __EXPORT int sbg_main(int argc, char *argv[])
{
	return SBG::main(argc, argv);
}
