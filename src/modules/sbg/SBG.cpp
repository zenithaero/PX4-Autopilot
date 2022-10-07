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
{
	// Reset timestamps
	_sbg_quat_data.timeStamp = 0;
	_sbg_nav_data.timeStamp = 0;
	_sbg_imu_data.timeStamp = 0;
	_sbg_gps_pos_data.timeStamp = 0;

	selector_status.timestamp = 0;
	sensor_selection.timestamp = 0;
	attitude.timestamp = 0;
	global_position.timestamp = 0;
	local_position.timestamp = 0;
	odometry.timestamp = 0;
}

// Static callback
SbgErrorCode SBG::onLogReceived(SbgEComHandle *pHandle, SbgEComClass msgClass, SbgEComMsgId msg, const SbgBinaryLogData *p_log_data, void *pUserArg)
{
	if (!SBG::instance) {
		return SBG_NULL_POINTER;
	}
	SBG::instance->onLog(msg, p_log_data);
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

	(void) _com_handle;
	(void) _interface;
	(void) _device_info;
	printf("[SBG CREATED]\n");
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
	errorCode = sbgEComInit(&_com_handle, &_interface);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to initialize the sbgECom library");
		return PX4_ERROR;
	}

	// Define callbacks for received data
	sbgEComSetReceiveLogCallback(&_com_handle, SBG::onLogReceived, NULL);

	// Configure output logs to 200 Hz
	errorCode = sbgEComCmdOutputSetConf(&_com_handle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_EKF_QUAT, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_EKF_QUAT");
		return PX4_ERROR;
	}

	errorCode = sbgEComCmdOutputSetConf(&_com_handle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_IMU_DATA, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_IMU_DATA");
		return PX4_ERROR;
	}

	errorCode = sbgEComCmdOutputSetConf(&_com_handle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_EKF_NAV, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_EKF_NAV");
		return PX4_ERROR;
	}

	errorCode = sbgEComCmdOutputSetConf(&_com_handle, SBG_ECOM_OUTPUT_PORT_A, SBG_ECOM_CLASS_LOG_ECOM_0, SBG_ECOM_LOG_GPS1_POS, SBG_ECOM_OUTPUT_MODE_MAIN_LOOP);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to configure output log SBG_ECOM_LOG_GPS1_POS");
		return PX4_ERROR;
	}

	// Display device information
	errorCode = sbgEComCmdGetInfo(&_com_handle, &_device_info);
	if(errorCode != SBG_NO_ERROR) {
		PX4_ERR("[SBG] Unable to get device information");
		return PX4_ERROR;
	}
	PX4_INFO("[SBG] Connected. SN: %.9u, version: %s\n", _device_info.serialNumber, SBG_E_COM_VERSION_STR);

	return PX4_OK;
}


void SBG::onLog(const SbgEComMsgId &msg, const SbgBinaryLogData *p_log_data)
{
	if (!p_log_data) {
		PX4_ERR("[SBG] Null log data exception");
		return;
	}
	bool pub_vehicle_attitude = false;
	bool pub_vehicle_local_position = false;
	bool pub_vehicle_global_position = false;
	bool pub_vehicle_odometry = false;

	switch (msg)
	{
	case SBG_ECOM_LOG_EKF_QUAT:
	{
		_sbg_quat_data = p_log_data->ekfQuatData;
		pub_vehicle_local_position = true;
		pub_vehicle_attitude = true;
		break;
	}
	case SBG_ECOM_LOG_IMU_DATA:
	{
		_sbg_imu_data = p_log_data->imuData;
		pub_vehicle_local_position = true;
		pub_vehicle_odometry = true;
		break;
	}
	case SBG_ECOM_LOG_EKF_NAV:
	{
		_sbg_nav_data = p_log_data->ekfNavData;
		pub_vehicle_local_position = true;
		pub_vehicle_global_position = true;
		pub_vehicle_odometry = true;
		break;
	}
	case SBG_ECOM_LOG_GPS1_POS:
	{
		_sbg_gps_pos_data = p_log_data->gpsPosData;
		pub_vehicle_local_position = true;
		pub_vehicle_global_position = true;
	}
	default:
		break;
	}

	if (pub_vehicle_attitude) publishVehicleAttitude();
	if (pub_vehicle_local_position) publishVehicleLocalPosition();
	if (pub_vehicle_global_position) publishVehicleGlobalPosition();
	if (pub_vehicle_odometry) publishVehicleOdometry();
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
	if (!_sbg_quat_data.timeStamp) return;

	attitude = {};
	attitude.q[0] = _sbg_quat_data.quaternion[0];
	attitude.q[1] = _sbg_quat_data.quaternion[1];
	attitude.q[2] = _sbg_quat_data.quaternion[2];
	attitude.q[3] = _sbg_quat_data.quaternion[3];
	//
	attitude.delta_q_reset[0] = 0;
	attitude.delta_q_reset[1] = 0;
	attitude.delta_q_reset[2] = 0;
	attitude.delta_q_reset[3] = 0;
	attitude.quat_reset_counter = 0;

	attitude.timestamp = hrt_absolute_time();
	_vehicle_attitude_pub.publish(attitude);

	publishVehicleOdometry();
}

void SBG::publishVehicleGlobalPosition()
{
	if (!_sbg_nav_data.timeStamp || !_sbg_gps_pos_data.timeStamp) return;

	global_position = {};
	global_position.lat = _sbg_nav_data.position[0];
	global_position.lon = _sbg_nav_data.position[1];
	global_position.alt = _sbg_nav_data.position[2];
	global_position.alt_ellipsoid = global_position.alt;

	global_position.delta_alt = 0;
	global_position.lat_lon_reset_counter = 0;
	global_position.alt_reset_counter = 0;

	global_position.eph = max(_sbg_gps_pos_data.latitudeAccuracy, _sbg_gps_pos_data.longitudeAccuracy);
	global_position.epv = _sbg_gps_pos_data.altitudeAccuracy;

	global_position.terrain_alt = 0;
	global_position.terrain_alt_valid = false;
	global_position.dead_reckoning = false;

	global_position.timestamp = hrt_absolute_time();
	global_position.timestamp_sample = hrt_absolute_time();


	// Init ref if needed
	if (!_pos_ref.isInitialized()) {
		setGlobalOrigin(global_position.lat, global_position.lon, global_position.alt);
	}

	_vehicle_global_position_pub.publish(global_position);
}

void SBG::publishVehicleLocalPosition()
{
	if (!_sbg_nav_data.timeStamp ||
	    !_sbg_gps_pos_data.timeStamp ||
	    !_sbg_imu_data.timeStamp ||
	    !_pos_ref.isInitialized()) return;

	local_position = {};

	local_position.xy_valid = true;
	local_position.z_valid = true;
	local_position.v_xy_valid = true;
	local_position.v_z_valid = true;

	// Position in local NED frame
	_pos_ref.project(
		_sbg_nav_data.position[0], _sbg_nav_data.position[1],
		local_position.x, local_position.y
	);
	local_position.z = _alt_ref - _sbg_nav_data.position[2];

	// Position reset delta
	local_position.delta_xy[0] = 0;
	local_position.delta_xy[1] = 0;
	local_position.xy_reset_counter = 0;

	local_position.delta_z = 0;
	local_position.z_reset_counter = 0;

	// Velocity in NED frame
	local_position.vx = _sbg_nav_data.velocity[0];
	local_position.vy = _sbg_nav_data.velocity[1];
	local_position.vz = _sbg_nav_data.velocity[2];
	local_position.z_deriv = local_position.vz;

	// Velocity reset delta
	local_position.delta_vxy[0] = 0;
	local_position.delta_vxy[1] = 0;
	local_position.vxy_reset_counter = 0;

	local_position.delta_vz = 0;

	// Acceleration in NED frame
	Quatf quat = Quatf(_sbg_quat_data.quaternion);
	Vector3f accel_body = Vector3f(
		_sbg_imu_data.accelerometers[0],
		_sbg_imu_data.accelerometers[1],
		_sbg_imu_data.accelerometers[2]);
	Vector3f accel_ned = quat.rotateVector(accel_body);

	local_position.ax = accel_ned(0);
	local_position.ay = accel_ned(1);
	local_position.az = accel_ned(2);

	Eulerf euler = Eulerf(quat);
	local_position.heading = euler(2);
	local_position.delta_heading = 0;
	local_position.heading_reset_counter = 0;
	local_position.heading_good_for_control = true;

	// Position of reference point (local NED frame origin) in global (GPS / WGS84) frame
	local_position.xy_global = true;
	local_position.z_global = true;
	local_position.ref_timestamp = _pos_ref.getProjectionReferenceTimestamp();
	local_position.ref_lat = _pos_ref.getProjectionReferenceLat();
	local_position.ref_lon = _pos_ref.getProjectionReferenceLon();
	local_position.ref_alt = _alt_ref;

	// Distance to surface
	local_position.dist_bottom = 0;
	local_position.dist_bottom_valid = false;
	local_position.dist_bottom_sensor_bitfield = vehicle_local_position_s::DIST_BOTTOM_SENSOR_NONE;

	local_position.eph = max(_sbg_gps_pos_data.latitudeAccuracy, _sbg_gps_pos_data.longitudeAccuracy);
	local_position.epv = _sbg_gps_pos_data.altitudeAccuracy;
	local_position.evh = 0;
	local_position.evv = 0;

	// Estimator specified vehicle limits
	local_position.vxy_max = 0;
	local_position.vz_max = 0;
	local_position.hagl_min = 0;
	local_position.hagl_max = 0;

	local_position.timestamp = hrt_absolute_time();
	local_position.timestamp_sample = hrt_absolute_time();
	_vehicle_local_position_pub.publish(local_position);

	publishVehicleOdometry();
}

void SBG::publishVehicleOdometry()
{
	if (!local_position.timestamp ||
            !attitude.timestamp ||
	    !_sbg_imu_data.timeStamp) return;

	odometry = {};

	odometry.local_frame = vehicle_odometry_s::LOCAL_FRAME_NED;
	odometry.x = local_position.x;
	odometry.y = local_position.y;
	odometry.z = local_position.z;

	odometry.q[0] = attitude.q[0];
	odometry.q[1] = attitude.q[1];
	odometry.q[2] = attitude.q[2];
	odometry.q[3] = attitude.q[3];

	odometry.q_offset[0] = 0;
	odometry.q_offset[1] = 0;
	odometry.q_offset[2] = 0;
	odometry.q_offset[3] = 0;

	memset(odometry.pose_covariance, 0, DIM(odometry.pose_covariance));

	odometry.velocity_frame = vehicle_odometry_s::LOCAL_FRAME_NED;

	odometry.vx = local_position.vx;
	odometry.vy = local_position.vy;
	odometry.vz = local_position.vz;

	odometry.rollspeed = _sbg_imu_data.gyroscopes[0];
	odometry.pitchspeed = _sbg_imu_data.gyroscopes[1];
	odometry.yawspeed = _sbg_imu_data.gyroscopes[2];

	memset(odometry.velocity_covariance, 0, DIM(odometry.velocity_covariance));

	odometry.reset_counter = 0;

	odometry.timestamp = hrt_absolute_time();
	odometry.timestamp_sample = hrt_absolute_time();
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
	SbgBinaryLogData p_log_data{};

	// Emulate attitude
	float w = 1.f; // rad/s
	float delta = 3.14 / 2;
	(void) delta;
	(void) ts;
	(void) w;
	// Eulerf euler = Eulerf(
	// 	delta * cosf(w * ts),
	// 	delta *cosf(w * ts),
	// 	delta *cosf(w * ts));
	Eulerf euler = Eulerf(
		0,
		0,
		0);
	Quatf quat = Quatf(euler);
	p_log_data.ekfQuatData.quaternion[0] = quat(0);
	p_log_data.ekfQuatData.quaternion[1] = quat(1);
	p_log_data.ekfQuatData.quaternion[2] = quat(2);
	p_log_data.ekfQuatData.quaternion[3] = quat(3);
	p_log_data.ekfQuatData.timeStamp = hrt_absolute_time();
	onLog(SBG_ECOM_LOG_EKF_QUAT, &p_log_data);

	// Emulate position
	double  lat = 48.8566;
	double lon = 2.3522;
	double altMSL = 120;

	p_log_data.ekfNavData.position[0] = lat;
	p_log_data.ekfNavData.position[1] = lon;
	p_log_data.ekfNavData.position[2] = altMSL + ts;

	if (_pos_ref.isInitialized()) {
		_pos_ref.reproject(
			0, ts,
			p_log_data.ekfNavData.position[0],
			p_log_data.ekfNavData.position[1]);
	}

	//
	p_log_data.ekfNavData.velocity[0] = 0.;
	p_log_data.ekfNavData.velocity[1] = 0.;
	p_log_data.ekfNavData.velocity[2] = 0.;
	p_log_data.ekfNavData.timeStamp = hrt_absolute_time();
	onLog(SBG_ECOM_LOG_EKF_NAV, &p_log_data);

	// Emulate gps pos
	p_log_data.gpsPosData.longitudeAccuracy = 1;
	p_log_data.gpsPosData.latitudeAccuracy = 1;
	p_log_data.gpsPosData.altitudeAccuracy = 1;
	p_log_data.gpsPosData.timeStamp = hrt_absolute_time();
	onLog(SBG_ECOM_LOG_GPS1_POS, &p_log_data);

	// Emulate acceleration
	p_log_data.imuData.accelerometers[0] = 0;
	p_log_data.imuData.accelerometers[1] = 0;
	p_log_data.imuData.accelerometers[2] = 9.81;
	p_log_data.imuData.gyroscopes[0] = 0;
	p_log_data.imuData.gyroscopes[1] = 0;
	p_log_data.imuData.gyroscopes[2] = 0;
	onLog(SBG_ECOM_LOG_IMU_DATA, &p_log_data);
}

void SBG::Run()
{
	if (should_exit()) {
		// _trigger_sub.unregisterCallback();
		exit_and_cleanup();
		return;
	}

	// Emulate data
	if (EMULATE) emulate();


	// Update global origin if needed
	updateGlobalOrigin();

	// Publish
	publishEstimatorSelectorStatus();
	publishSensorSelection();
	publishVehicleAttitude();
	publishVehicleLocalPosition();
	publishVehicleGlobalPosition();
	publishVehicleOdometry();
	publishWindEstimate();

	// Temp
	printf("[SBG loop]\n");
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
SBG IMU Driver
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
