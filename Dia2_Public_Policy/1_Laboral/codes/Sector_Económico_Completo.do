** Sectores Completos a CIIU Rev.4
* Actividad económica
gen sector1r4=.
replace sector1r4=1 if p506r4>=100 & p506r4<300
replace sector1r4=2 if p506r4>=300 & p506r4<500
replace sector1r4=3 if p506r4>=500 & p506r4<1000
replace sector1r4=4 if p506r4>=1000 & p506r4<3500
replace sector1r4=5 if p506r4>=3500 & p506r4<3600
replace sector1r4=6 if p506r4>=3600 & p506r4<4100
replace sector1r4=7 if p506r4>=4100 & p506r4<4500
replace sector1r4=8 if p506r4>=4500 & p506r4<4900
replace sector1r4=9 if p506r4>=4900 & p506r4<5500
replace sector1r4=10 if p506r4>=5500 & p506r4<5800
replace sector1r4=11 if p506r4>=5800 & p506r4<6400
replace sector1r4=12 if p506r4>=6400 & p506r4<6800
replace sector1r4=13 if p506r4>=6800 & p506r4<6900
replace sector1r4=14 if p506r4>=6900 & p506r4<7700
replace sector1r4=15 if p506r4>=7700 & p506r4<8400
replace sector1r4=16 if p506r4>=8400 & p506r4<8500
replace sector1r4=17 if p506r4>=8500 & p506r4<8600
replace sector1r4=18 if p506r4>=8600 & p506r4<9000
replace sector1r4=19 if p506r4>=9000 & p506r4<9400
replace sector1r4=20 if p506r4>=9400 & p506r4<9700
replace sector1r4=21 if p506r4>=9700 & p506r4<9900
replace sector1r4=22 if p506r4==9900 

#delimit ;
lab def sector1r4 
1 "Agricultura, ganadería y silvicultura" 2 "Pesca y acuicultura" 3 " Explotación de minas y canteras" 4 "Industrias manufactureras" 5 "Suministro de electricidad, gas, vapor y aire acondicionado"
6 "Suministro de agua; evacuación de aguas residuales, gestión de desechos y descontaminación" 7 "Construcción" 
8 "Comercio al por mayor y al por menor; reparación de vehículos automotores y motocicletas" 9 "Transporte y almacenamiento" 
10 "Actividades de alojamiento y de servicio de comidas" 11 "Información y comunicaciones" 12 "Actividades financieras y de seguros" 13 " Actividades inmobiliarias"
14 "Actividades profesionales, científicas y técnicas" 15 "Actividades de servicios administrativos y de apoyo" 16 "Administración pública y defensa; planes de seguridad social de afiliación obligatoria " 
17 "Enseñanza" 18 "Actividades de atención de la salud humana y de asistencia social" 19 " Actividades artísticas, de entretenimiento y recreativas" 20 "Otras actividades de servicios"
21 "Actividades de los hogares como empleadores; actividades no diferenciadas de los hogares como productores de bienes y servicios para uso propio"
22 "Actividades de organizaciones y órganos extraterritoriales"; 
#delimit cr
lab val sector1r4 sector1r4 

* Sector económico
gen sector2r4_2=.
replace sector2r4_2=1 if p506r4>=100 & p506r4<300
replace sector2r4_2=2 if p506r4>=300 & p506r4<500
replace sector2r4_2=3 if p506r4>=500 & p506r4<1000
replace sector2r4_2=4 if p506r4>=1000 & p506r4<3500
replace sector2r4_2=5 if p506r4>=3500 & p506r4<4100
replace sector2r4_2=6 if p506r4>=4100 & p506r4<4500
replace sector2r4_2=7 if p506r4>=4500 & p506r4<4900
replace sector2r4_2=8 if p506r4>=4900 & p506r4<5500
replace sector2r4_2=9 if p506r4>=5500 & p506r4<9998
#delimit ;
lab def sector2r4_2 
1 "Agricultura, ganadería y silvicultura" 
2 "Pesca y acuicultura" 3 " Explotación de minas y canteras" 4 "Industrias manufactureras" 
5 "Suministro de electricidad, gas y agua" 6 "Construcción" 7 "Comercio" 8 "Transporte y almacenamiento" 9 "Servicios";
#delimit cr
lab val sector2r4_2 sector2r4_2

* Sector económico agregado
gen sector3r4=. 
replace sector3r4=1 if sector2r4_2==1
replace sector3r4=2 if sector2r4_2==2
replace sector3r4=3 if sector2r4_2==3
replace sector3r4=4 if sector2r4_2==4
replace sector3r4=5 if sector2r4_2==5 | sector2r4_2==8 | sector2r4_2==9
replace sector3r4=6 if sector2r4_2==6
replace sector3r4=7 if sector2r4_2==7
#delimit ;
lab def sector3r4 1 "Agricultura" 2 "Pesca y acuicultura" 3 " Minería" 4 "Manufactura" 
5 "Servicios" 6 "Construcción" 7 "Comercio";
#delimit cr
lab val sector3r4 sector3r4








