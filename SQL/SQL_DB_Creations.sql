use cas_delay_test;
create table state(
	lat double,
    lon double,
    Heading double,
	V_meters double,
	V_knots double ,
	timeStamp double ,
	ACID varchar(100),
	ExpNum int,
	Category varchar(100),
	ExpType varchar(200),
	ax_spacing double
);

create table note_events(
	event varchar(500),
	timeStamp double,
	ACID varchar(100),
	ExpNum int,
	Category varchar(100),
	ExpType varchar(200),
	ax_spacing double
);

create table ev_specific(
	TOI double,
	TOI_Dist double,
	timeStamp double,
	ax_spacing double,
	ACID varchar(100),
	ExpNum int,
	Category varchar(100),
	ExpType varchar(200)
);

create table dubins(
path_type varchar(100),
bez_intx double,
bez_inty double,
nom_intx double,
nom_inty double,
bez_t double,
intersect_heading double,
h double,
k double,
tr double,
bank_angle double,
velocity double,
timeStamp double,
ACID varchar(100),
ExpNum int,
Category varchar(100),
ExpType varchar(200),
ax_spacing double
);

create table bez(
	p0x double,
	p0y double,
	p1x double,
	p1y double,
	p2x double,
	p2y double,
	length double,
	Bez_ID varchar(100),
	TOA double,
	orig_TOA double,   
	travel_time double,    
	timeStamp double,     
	velocity double,      
	ACID varchar(100),
	ExpNum int,
	Category varchar(100),
	ExpType varchar(200),
	ax_spacing double
);

create table aircraft(
ACID varchar(100),
ExpNum int,
Category varchar(100),
dt double,
ExpType varchar(200),
ax_spacing double
);
create table delay(
delay double,
nom_time double,
am_time double,
point_dist double,
start_lat double,
start_lon double,
end_lat double,
end_lon double,
ACID varchar(100),
ExpNum int,
Category varchar(100),
ExpType varchar(200),
ax_spacing double

);