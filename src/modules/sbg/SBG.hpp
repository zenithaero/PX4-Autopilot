#pragma once

#include <lib/mathlib/mathlib.h>
#include <lib/parameters/param.h>
#include <px4_platform_common/px4_config.h>
#include <px4_platform_common/defines.h>
#include <px4_platform_common/module.h>
#include <px4_platform_common/module_params.h>
#include <px4_platform_common/posix.h>
#include <uORB/Publication.hpp>
#include <uORB/Subscription.hpp>
#include <uORB/SubscriptionCallback.hpp>

#include <px4_platform_common/px4_work_queue/ScheduledWorkItem.hpp>
#include <lib/geo/geo.h>
#include <sbgEComLib.h>
// Subscriptions
#include <uORB/topics/vehicle_command.h>
#include <uORB/topics/estimator_status.h>
// Raw publications
#include <uORB/topics/sbg_ekf_nav.h>
#include <uORB/topics/sbg_ekf_quat.h>
// Estimator publications
#include <uORB/topics/estimator_selector_status.h>
#include <uORB/topics/estimator_status.h>
#include <uORB/topics/parameter_update.h>
#include <uORB/topics/sensor_selection.h>
#include <uORB/topics/sensors_status_imu.h>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_local_position.h>
#include <uORB/topics/vehicle_global_position.h>
#include <uORB/topics/vehicle_odometry.h>
#include <uORB/topics/wind.h>


#define DIM(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

class SBG : public ModuleBase<SBG>, public ModuleParams, public px4::ScheduledWorkItem
{
public:
	static SBG *instance;

	SBG();
	~SBG() override = default;

	/** @see ModuleBase */
	static int task_spawn(int argc, char *argv[]);

	/** @see ModuleBase */
	static int custom_command(int argc, char *argv[]);

	/** @see ModuleBase */
	static int print_usage(const char *reason = nullptr);

	bool init();

private:
	// Sbg interface
	SbgEComHandle _comHandle;
	SbgInterface _interface;
	SbgEComDeviceInfo _deviceInfo;

	// Previous messages
	estimator_selector_status_s selector_status{};
	sensor_selection_s sensor_selection{};
	vehicle_attitude_s attitude{};
	vehicle_global_position_s global_position{};
	vehicle_local_position_s local_position{};
	vehicle_odometry_s odometry{};
	wind_s wind{};


	void Run() override;


	// Subscriptions
	uORB::Subscription _vehicle_command_sub{ORB_ID(vehicle_command)};
	uORB::Subscription _estimator_status_sub{ORB_ID(estimator_status), 0};

	// Raw publications
	uORB::Publication<sbg_ekf_nav_s>	_sbg_ekf_nav_pub{ORB_ID(sbg_ekf_nav)};
	uORB::Publication<sbg_ekf_quat_s>	_sbg_ekf_quat_pub{ORB_ID(sbg_ekf_quat)};

	// Estimator publications
	uORB::Publication<estimator_selector_status_s> _estimator_selector_status_pub{ORB_ID(estimator_selector_status)};
	uORB::Publication<sensor_selection_s>          _sensor_selection_pub{ORB_ID(sensor_selection)};
	uORB::Publication<vehicle_attitude_s>          _vehicle_attitude_pub{ORB_ID(vehicle_attitude)};
	uORB::Publication<vehicle_global_position_s>   _vehicle_global_position_pub{ORB_ID(vehicle_global_position)};
	uORB::Publication<vehicle_local_position_s>    _vehicle_local_position_pub{ORB_ID(vehicle_local_position)};
	uORB::Publication<vehicle_odometry_s>          _vehicle_odometry_pub{ORB_ID(vehicle_odometry)};
	uORB::Publication<wind_s>             _wind_pub{ORB_ID(wind)};

	// Data
	MapProjection _pos_ref{}; // WGS-84 position latitude and longitude of the origin
	float _alt_ref{0.0f}; // WGS-84 height (m)
	SbgLogEkfQuatData _sbgQuatData;
	SbgLogEkfNavData _sbgNavData;
	SbgLogImuData _sbgImuData;

	// Serial logs callback
	static SbgErrorCode onLogReceived(SbgEComHandle *pHandle, SbgEComClass msgClass, SbgEComMsgId msg, const SbgBinaryLogData *pLogData, void *pUserArg);


	int connect();
	void onLog(const SbgEComMsgId &msg, const SbgBinaryLogData *pLogData);

	void emulate();
	void updateGlobalOrigin();
	void setGlobalOrigin(const double latitude, const double longitude, const float altitude);
	void getGlobalOrigin(double &latitude, double &longitude, float &origin_alt) const;

	void publishEstimatorSelectorStatus();
	void publishSensorSelection();
	void publishVehicleAttitude();
	void publishVehicleLocalPosition();
	void publishVehicleGlobalPosition();
	void publishVehicleOdometry();
	void publishWindEstimate();

	// uORB::SubscriptionCallbackWorkItem _trigger_sub{this, ORB_ID(camera_trigger)};
	// param_t _p_cam_cap_fback;
	// int32_t _cam_cap_fback{0};
};
