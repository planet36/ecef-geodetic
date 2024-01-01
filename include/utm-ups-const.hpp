// SPDX-FileCopyrightText: Steven Ward
// SPDX-License-Identifier: OSL-3.0

/// Transverse Mercator conversions
/**
\file
\author Steven Ward
\sa https://geographiclib.sourceforge.io/C++/doc/transversemercator.html#tmseries
\sa https://geographiclib.sourceforge.io/html/tmseries30.html
*/

// Generated by "generate-series-approx-coeff-tm.cal" on (@1695775826) Tue Sep 26 20:50:26 2023

#pragma once

#include <concepts>
#include <valarray>

/// Geodetic constants
namespace Geodetic_const
{
inline constexpr int min_lat = -90; // (degrees)
inline constexpr int max_lat =  90; // (degrees)

inline constexpr int min_lon = -180; // (degrees)
inline constexpr int max_lon =  180; // (degrees)
}

/// transverse Mercator constants
namespace TM_const
{
template <std::floating_point T>
inline constexpr T A = T{6.367449145823415309285117466604332246914402102957570737850939755154909040e6L};

template <std::floating_point T>
const std::valarray<T> alpha{
	T{8.377318206244698303199297980264726661812122018748626612607773594358891396e-4L},
	T{7.608527773572489156370563711094500187721012371404860521900273049502074293e-7L},
	T{1.197645503242492100577207243355655395550437605188580688459348441237652293e-9L},
	T{2.429170680397313150152152022748612877329876862299544259783100295210212997e-12L},
	T{5.711818369154105277380619294463053128246705077533709943589895010989052231e-15L},
	T{1.479998027052620835427117989571514820485499047805091801402854352965848021e-17L},
	T{4.107687520524031171552835537199016012142739701011294776412550867990043025e-20L},
	T{1.199991067742438187177239677609345487818171856341259153750529821495206319e-22L},
	T{3.647332481594558850072446840415251984173752875298642123816051858107440385e-25L},
	T{1.144173305631218813924026394149837692176530386128546857016161361233191053e-27L},
	T{3.682951266782108579534818990177535732300575780601610932702049139434209273e-30L},
	T{1.211146842124174288757500192264621126482155774456938094069349746402656385e-32L},
	T{4.055501161666398024114556580356853347077563706087579913014396199218313161e-35L},
	T{1.379124152150737733828358112131123588786699856847831539980836602969844881e-37L},
	T{4.752991037031996160205450888568786453115865636538277055628433708396095216e-40L},
	T{1.657298494311255715533736831465587148652815710563912816177206627303954602e-42L},
	T{5.838493517162075915321038121586944604957851178340818171000806164916457555e-45L},
	T{2.075701948517379047265234336864650863166157124729579301294062547513884232e-47L},
	T{7.439969154510613375982513107934654186901610689856023478526392553971643537e-50L},
	T{2.686339863387493904934826507092169663426188251516305473720982291729650018e-52L},
	T{9.764001230453233942575134002153376866486789901167781481457636961986479530e-55L},
	T{3.570321206051312165054476853094037086250772240100290340633488621872239122e-57L},
	T{1.312712706842644602552675299165171972984485576828699995115322328956856171e-59L},
	T{4.850810596036506891430184403555079560206020505312421709357098724936752462e-62L},
	T{1.800796985553215411797627913795817543416022954448591925059689104071487968e-64L},
	T{6.713749096025871466549483854195152406665960001808668361503829634658086007e-67L},
	T{2.512910071379620997734287466707960938776293687626807666774328519878020824e-69L},
	T{9.440099084830310047202484218470808788112057417924647790521459862899427575e-72L},
	T{3.557861351934538984949889172962064409819792073546146677784488239810091291e-74L},
	T{1.376701097054089722281774223709505617466106944371738402676856193827844269e-76L},
};

template <std::floating_point T>
const std::valarray<T> beta{
	T{8.377321640579486770657338660487753534316992195008696388559132269213256879e-4L},
	T{5.905870152220365181507106735379055148798627384800312313658117065473551841e-8L},
	T{1.673482665343824927368284807979927027879661080730636441133772729883005373e-10L},
	T{2.164798110490386180776881159858920963592981699109268506764648203859954329e-13L},
	T{3.787930968839601126376873899353149941085908277282374724195568410190961853e-16L},
	T{7.236769287965757922603588687959988305076578982215997474249259064535054420e-19L},
	T{1.493454494878010126640350065464049110696405049691319428631175307373309155e-21L},
	T{3.253843083999393961181907869952142417582220825644524824173728914231013795e-24L},
	T{7.391210921735273959428771950905376597587691323070994340634200582795664753e-27L},
	T{1.735023326722465427935103103705924909828754571441064206776786823766958485e-29L},
	T{4.182537520202514777241166498102002896477051033005346930697132933238329272e-32L},
	T{1.030635528246784830298318470313701687367395044892107846282983532932278071e-34L},
	T{2.586855862336712844488060328551228113531925326246712636720834176042765767e-37L},
	T{6.595581491129959501196193345536390786537947617002899231944391450907125082e-40L},
	T{1.704535565040829396197199359608901393209230368975798087599887576323716786e-42L},
	T{4.457303859983084977645606816412090694161752476824653171653148779445781392e-45L},
	T{1.177690003743177367822788449260089778381587184591130224316104012470371416e-47L},
	T{3.140279406512039704932874366065560229570835272163704890335766550496281182e-50L},
	T{8.442161100395104591035460015605118765723974102803354284797374508514832345e-53L},
	T{2.286245418731765634549755039310382224576603547042881344705246688564707250e-55L},
	T{6.232540747460568523109358170924943102454285947594407616949627333415884021e-58L},
	T{1.709279527843028084943935714314144514651053615990531651016978042248800282e-60L},
	T{4.713408208684100282867203266395590692431362749447027795703736066782801198e-63L},
	T{1.306257911556640840711958626730042557734336769388145526345388350088270817e-65L},
	T{3.636783238968222285218207219725748311174559182549775257549135788249988619e-68L},
	T{1.016819966035676279114678511754644843867486828625291337952663815948065534e-70L},
	T{2.854106423931739540306310710097106881690113131647859332519048158865055139e-73L},
	T{8.040296766555112821029670182930508298889196186739844989495843298065690370e-76L},
	T{2.272830222485031211789453855604682309392148483051170710498105085148619356e-78L},
	T{6.457527744777211098050789700274549842364047682210954928686657110916244540e-81L},
};

}

/// polar stereographic constants
namespace PS_const
{
/// (1-f)*exp(e*atanh(e)) == sqrt((1+e)^(1+e) * (1-e)^(1-e))
template <std::floating_point T>
inline constexpr T k90 = T{1.003356555249315433614531911263106313397037011923772178096177286730040894L};
}

/// Universal Transverse Mercator constants
namespace UTM_const
{
inline constexpr int min_zone = 1;
inline constexpr int max_zone = 60;
inline constexpr int zone_width = 6; // (degrees)
inline constexpr int half_zone_width = zone_width / 2; // (degrees)

/// central meridian scale factor (unitless)
template <std::floating_point T>
inline constexpr T k0 = T{0.9996L};

inline constexpr int false_easting = 500'000; // easting at central meridian (meters)
inline constexpr int min_easting = 100'000; // (meters)
inline constexpr int max_easting = 900'000; // (meters)

inline constexpr int false_northing[2] = {10'000'000, 0}; // [south, north] northing at equator (meters)
inline constexpr int min_northing[2] = {300'000, 0}; // [south, north] (meters)
inline constexpr int max_northing[2] = {10'000'000, 9'700'000}; // [south, north] (meters)
}

/// Universal Polar Stereographic constants
namespace UPS_const
{
/// central meridian scale factor (unitless)
template <std::floating_point T>
inline constexpr T k0 = T{0.994L};

// Note: The easting and northing are symmetric and have the same constraints.
inline constexpr int false_easting_northing = 2'000'000; // easting/northing at pole (meters)
inline constexpr int min_easting_northing[2] = {800'000, 1'300'000}; // [south, north] (meters)
inline constexpr int max_easting_northing[2] = {3'200'000, 2'700'000}; // [south, north] (meters)
}

/// Military Grid Reference System constants
namespace MGRS_const
{
inline constexpr int num_lat_bands = 20;

inline constexpr int min_lat = -80; // (degrees)
inline constexpr int max_lat = 84; // (degrees)
inline constexpr int lat_band_height = 8; // (degrees)

inline constexpr int meters_per_square = 100'000; // (meters/square)
inline constexpr int northing_period_squares = 20; // (squares)
inline constexpr int northing_period_meters = meters_per_square * northing_period_squares; // (meters)
}
