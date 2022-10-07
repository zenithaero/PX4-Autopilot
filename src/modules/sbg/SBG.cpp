#include "SBG.hpp"
// #include <math.h>
#include <matrix/math.hpp>

using namespace time_literals;
using namespace matrix;

#define EKF2_MAX_INSTANCES 9
#define EMULATE true

SBG* SBG::instance = nullptr;

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
	int ret = EMULATE ? 0 : connect();
	if (ret) {
		PX4_ERR("[SBG] Connection error");
		return false;
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
	ScheduleOnInterval(100_ms);


	return true;
}

int SBG::connect()
{
	SbgErrorCode errorCode;
	errorCode = sbgInterfaceSerialCreate(&_interface, "/dev/tty0", 9600);
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
	if (!pLogData) {
		PX4_ERR("[SBG] Null log data exception");
		return;
	}
	bool pubVehicleAttitude = false;
	bool pubVehicleLocalPosition = false;
	bool pubVehicleGlobalPosition = false;
	bool pubVehicleOdometry = false;

	switch (msg)
	{
	case SBG_ECOM_LOG_EKF_QUAT:
	{
		_sbgQuatData = pLogData->ekfQuatData;
		pubVehicleAttitude = true;
		break;
	}
	case SBG_ECOM_LOG_IMU_DATA:
	{
		_sbgImuData = pLogData->imuData;
		pubVehicleLocalPosition = true;
		pubVehicleOdometry = true;
		break;
	}
	case SBG_ECOM_LOG_EKF_NAV:
	{
		_sbgNavData = pLogData->ekfNavData;
		pubVehicleLocalPosition = true;
		pubVehicleGlobalPosition = true;
		pubVehicleOdometry = true;
		break;
	}
	default:
		break;
	}

	if (pubVehicleAttitude) publishVehicleAttitude();
	if (pubVehicleLocalPosition) publishVehicleLocalPosition();
	if (pubVehicleGlobalPosition) publishVehicleGlobalPosition();
	if (pubVehicleOdometry) publishVehicleOdometry();
}

void SBG::publishEstimatorSelectorStatus()
{
	// Retrieve ids
	estimator_status_s estimator_status{};
	_estimator_status_sub.copy(&estimator_status);

	selector_status = {};
	selector_status.primary_instance = 0;
	selector_status.instances_available = 1;
	selector_status.instance_changed_count = 0;
	selector_status.last_instance_change = 0;
	selector_status.accel_device_id = estimator_status.accel_device_id;
	selector_status.baro_device_id = estimator_status.baro_device_id;
	selector_status.gyro_device_id = estimator_status.gyro_device_id;
	selector_status.mag_device_id = estimator_status.mag_device_id;
	selector_status.gyro_fault_detected = 0;
	selector_status.accel_fault_detected = 0;

	for (size_t i = 0; i < DIM(selector_status.healthy); i++) {
		selector_status.combined_test_ratio[i] = 0;
		selector_status.relative_test_ratio[i] = 0;
		selector_status.healthy[i] = true;
	}

	for (size_t i = 0; i < DIM(selector_status.accumulated_gyro_error); i++) {
		selector_status.accumulated_gyro_error[i] = 0;
		selector_status.accumulated_accel_error[i] = 0;
	}

	selector_status.timestamp = hrt_absolute_time();
	_estimator_selector_status_pub.publish(selector_status);
}

void SBG::publishSensorSelection()
{
	// Retrieve ids
	estimator_status_s estimator_status{};
	_estimator_status_sub.copy(&estimator_status);

	sensor_selection = {};
	sensor_selection.accel_device_id = estimator_status.accel_device_id;
	sensor_selection.gyro_device_id = estimator_status.gyro_device_id;
	sensor_selection.timestamp = hrt_absolute_time();
	_sensor_selection_pub.publish(sensor_selection);
}

void SBG::publishVehicleAttitude()
{
	attitude = {};
	attitude.q[0] = _sbgQuatData.quaternion[0];
	attitude.q[1] = _sbgQuatData.quaternion[1];
	attitude.q[2] = _sbgQuatData.quaternion[2];
	attitude.q[3] = _sbgQuatData.quaternion[3];
	//
	attitude.delta_q_reset[0] = 0;
	attitude.delta_q_reset[1] = 0;
	attitude.delta_q_reset[2] = 0;
	attitude.delta_q_reset[3] = 0;
	attitude.quat_reset_counter = 0;

	attitude.timestamp = hrt_absolute_time();
	_vehicle_attitude_pub.publish(attitude);
}

void SBG::publishVehicleLocalPosition()
{
	local_position = {};

	local_position.timestamp = hrt_absolute_time();
	_vehicle_local_position_pub.publish(local_position);
}

void SBG::publishVehicleGlobalPosition()
{
	global_position = {};

	global_position.lat = _sbgNavData.position[0];
	global_position.lon = _sbgNavData.position[1];
	global_position.alt = _sbgNavData.position[2];
	global_position.alt_ellipsoid = global_position.alt;

	global_position.delta_alt = 0;
	global_position.lat_lon_reset_counter = 0;
	global_position.alt_reset_counter = 0;

	global_position.eph = 0;
	global_position.epv = 0;

	global_position.terrain_alt = 0;
	global_position.terrain_alt_valid = false;
	global_position.dead_reckoning = false;

	global_position.timestamp = hrt_absolute_time();
	global_position.timestamp_sample = hrt_absolute_time();
	_vehicle_global_position_pub.publish(global_position);
}

void SBG::publishVehicleOdometry()
{
	odometry = {};

	odometry.timestamp = hrt_absolute_time();
	_vehicle_odometry_pub.publish(odometry);
}


void SBG::publishWindEstimate()
{
	wind = {};
	wind.windspeed_north = 0;
	wind.windspeed_east = 0;
	wind.variance_north = 0;
	wind.variance_east = 0;
	wind.tas_innov = 0;
	wind.tas_innov_var = 0;
	wind.beta_innov = 0;
	wind.beta_innov_var = 0;

	wind.timestamp_sample = hrt_absolute_time();
	wind.timestamp = hrt_absolute_time();
	_wind_pub.publish(wind);
}


void SBG::setGlobalOrigin(const double latitude, const double longitude, const float altitude)
{
	_pos_ref.initReference(latitude, longitude);
	_alt_ref = altitude;
}

void SBG::getGlobalOrigin(double &latitude, double &longitude, float &origin_alt) const
{
	latitude = _pos_ref.getProjectionReferenceLat();
	longitude = _pos_ref.getProjectionReferenceLon();
	origin_alt  = _alt_ref;
}

void SBG::updateGlobalOrigin()
{
	if (!_vehicle_command_sub.updated()) return;
	vehicle_command_s vehicle_command;
	if (!_vehicle_command_sub.update(&vehicle_command)) return;
	if (vehicle_command.command == vehicle_command_s::VEHICLE_CMD_SET_GPS_GLOBAL_ORIGIN) {
		double latitude = vehicle_command.param5;
		double longitude = vehicle_command.param6;
		float altitude = vehicle_command.param7;
		setGlobalOrigin(latitude, longitude, altitude);
		getGlobalOrigin(latitude, longitude, altitude);
		PX4_INFO("New NED origin (LLA): %3.10f, %3.10f, %4.3f\n", latitude, longitude, static_cast<double>(altitude));
	}
}

void SBG::emulate()
{
	float ts = hrt_absolute_time() / 1e6f;
	SbgBinaryLogData pLogData{};

	// Emulate attitude
	float w = 1.f; // rad/s
	float delta = 3.14 / 2;
	Eulerf euler = Eulerf(
		delta * cosf(w * ts),
		delta *cosf(w * ts),
		delta *cosf(w * ts));
	Quatf quat = Quatf(euler);
	pLogData.ekfQuatData.quaternion[0] = quat(0);
	pLogData.ekfQuatData.quaternion[1] = quat(1);
	pLogData.ekfQuatData.quaternion[2] = quat(2);
	pLogData.ekfQuatData.quaternion[3] = quat(3);
	onLog(SBG_ECOM_LOG_EKF_QUAT, &pLogData);

	// SBG_ECOM_LOG_EKF_QUAT
	// SBG_ECOM_LOG_IMU_DATA
	// SBG_ECOM_LOG_EKF_NAV

	// Emulate position
	double  lat = 48.8566;
	double lon = 2.3522;
	double altMSL = 120;
	//
	pLogData.ekfNavData.position[0] = lat;
	pLogData.ekfNavData.position[1] = lon;
	pLogData.ekfNavData.position[2] = altMSL;
	//
	pLogData.ekfNavData.velocity[0] = 0.;
	pLogData.ekfNavData.velocity[1] = 0.;
	pLogData.ekfNavData.velocity[2] = 0.;
	onLog(SBG_ECOM_LOG_EKF_NAV, &pLogData);
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

	// Emulate data
	if (EMULATE) emulate();




	// Publish
	publishEstimatorSelectorStatus();
	publishSensorSelection();
	publishVehicleAttitude();
	publishVehicleLocalPosition();
	publishVehicleGlobalPosition();
	publishVehicleOdometry();
	publishWindEstimate();
	printf("[SBG loop]\n");

	// Publish estimator status
	// estimator_selector_status_s selector_status{};
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
