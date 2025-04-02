
/*==============================================================================
Proyecto:				Mercado Laboral Peruano
Autores:				Mario Huarancca Bellido	(mario.huarancca@bcrp.gob.pe)
Institución: 			BCRP
----------------------------------------------------------------------
Fecha de creación: 		05/03/2025
Fecha de modificación:	06/03/2025
Fuente:					INEI - ENAHO 2017 y 2022
Producto:				Tablas y Gráficos
==============================================================================*/

clear all
global data 	"C:\Users\Curso\Documents\Politicas publicas\Dia2_Public_Policy\1_Laboral\data"
global codes 	"C:\Users\Curso\Documents\Politicas publicas\Dia2_Public_Policy\1_Laboral\codes"
global output 	"C:\Users\Curso\Documents\Politicas publicas\Dia2_Public_Policy\1_Laboral\output"

**# Estructura del Mercado Laboral Peruano

	*** Unión de BDs
	use "${data}/enaho01a-2017-500", clear
	keep conglome vivienda hogar codperso ubigeo dominio estrato p204 p205 p206 p207 p501 p507 p509 p510 i524a1 p506r4 ///
		 d529t i530a d536 i538a1 d540t i541a d543 d544t ocu500 fac500a
	gen year = 2017
	append using "${data}/enaho01a-2022-500", keep(conglome vivienda hogar codperso ubigeo dominio estrato p204 p205 p206 p207 p501 p507 p509 p510 i524a1 p506r4 d529t ///
												   i530a d536 i538a1 d540t i541a d543 d544t ocu500 fac500a)
	replace year = 2022 if year==.

	*** Filtro del residente habitual
	keep if ((p204==1 & p205==2) | (p204==2 & p206==1)) & p501!=.

	*** Variables

		** Población en Edad de Trabajar (PET)
		gen pet = 1

		** Población Económicament e Activa (PEA)
		gen pea = (ocu500<=2)

		** PEA Ocupada
		gen peao = (ocu500==1)

		** PEA Desempleada
		gen pead = (ocu500==2)

		** Población Económicamente Inactiva (PEI)
		gen pei = (ocu500==3 | ocu500==4)

		** Tasa de ocupación
		gen tasa_ocup = 0
		replace tasa_ocup = 1 if ocu500==1

		** Sexo
		gen sexo = (p207==1)
		label define sexo 0 "Mujer" 1 "Hombre"
		label values sexo sexo

		** Área geográfica
		gen area = (estrato<=5)
		label define area 0 "Rural" 1 "Urbano"
		label values area area

		** Categoría ocupacional
		gen cat_ocup = .
		replace cat_ocup = 1 if (p507== 1 | p509==1)
		replace cat_ocup = 2 if (p507==3 & (p510>=3 & p510<=7) )
		replace cat_ocup = 3 if (p507==3 & (p510==1 | p510==2))
		replace cat_ocup = 4 if ((p507==4) & (p510==3 | p510==4 | p510==5 | p510==6 | p510==7)) 
		replace cat_ocup = 5 if (p507==4 & (p510==1 | p510==2) )
		replace cat_ocup = 6 if (p507==2 | (p507==1 & p509==2))
		replace cat_ocup = 7 if (p507==5 | p507==7)
		replace cat_ocup = 8 if (p507==6)
		label define cat_ocup 	1 "Empleador" 2 "Empleado privado" 3 "Empleado público" 4 "Obrero privado" 5 "Obrero público" ///
								6 "Independiente" 7 "Trab. familiar no remunerado" 8 "Trabajador del hogar"
		label values cat_ocup cat_ocup
		
		** Actividad económica
		do "${codes}/Sector_Económico_Completo.do"

		** Ingreso laboral
		egen ingreso = rsum(i524a1 d529t i530a d536 i538a1 d540t i541a d543 d544t) if ocu500==1
		replace ingreso = ingreso/12

	*** Tabulados

		** Cuadro 2.1 - Distribución de la PET según condición de actividad
		table () (year) [iw=fac500a], stat(sum pet pea peao pead pei) nformat(%10.0f) nototal

		** Gráfico 2.3 - Tasa de ocupación según características demográficas
		table (sexo) (year) [iw=fac500a], stat(mean tasa_ocup) nformat(%10.3f) total(year)
		table (area) (year) [iw=fac500a], stat(mean tasa_ocup) nformat(%10.3f) total(year)

		** Gráfico 2.10 - PEA ocupada según categoría ocupacional
		table (cat_ocup) (year) if ocu500==1 [iw=fac500a], stat(percent, across(cat_ocup)) nformat(%3.1f)

		** Salario promedio
		
			* Categoría ocupacional
			table (cat_ocup) (year) if ocu500==1 & (ingreso>0 & !missing(ingreso)) [iw=fac500a], stat(mean ingreso) total(year) nformat(%10.0f)
			collect export "${output}/tables/Cuadros.xlsx", sheet("C1", replace) cell(B4) modify 
		** El missing es necesario para no contar los que no registran

			* Actividad económica
			table (sector3r4) (year) if ocu500==1 & (ingreso>0 & !missing(ingreso)) [iw=fac500a], stat(mean ingreso) total(year) nformat(%10.0f)

**# Remuneración Mínima Vital (RMV)

	use "${data}/Sample_Empleo_junio_2023.dta", clear

	*** Variables

		** Tamaño de empresas
		gen tamaño = .
		replace tamaño = 1 if tamaão_empresa=="DE 1 A 10 TRABAJADORES"
		replace tamaño = 2 if tamaão_empresa=="DE 11 A 50 TRABAJADORES" | tamaão_empresa=="DE 51 A 100 TRABAJADORES"
		replace tamaño = 3 if tamaão_empresa=="DE 101 A 500 TRABAJADORES" | tamaão_empresa=="DE 501 A MÃS TRABAJADORES"
		label variable tamaño "Tamaño de empresas"
		label define tamaño 1 "Microempresa (1 a 10 trabajadores)" 2 "Pequeña empresa (11 a 100 trabajadores)" ///
							3 "Mediana y Gran empresa (más de 100 trabajadores)"
		label values tamaño tamaño

		** Sectores económicos desagregado
		gen sector = .
		replace sector = 1 if actividad_economica=="AGRICULTURA, GANADERÃA, CAZA Y SILVICULTURA"
		replace sector = 2 if actividad_economica=="PESCA"
		replace sector = 3 if actividad_economica=="EXPLOTACIÃN DE MINAS Y CANTERAS"
		replace sector = 4 if actividad_economica=="INDUSTRIAS MANUFACTURERAS"
		replace sector = 5 if actividad_economica=="SUMINISTRO DE ELECTRICIDAD, GAS Y AGUA"
		replace sector = 6 if actividad_economica=="CONSTRUCCIÃN"
		replace sector = 7 if actividad_economica=="COMERCIO AL POR MAYOR Y AL POR MENOR, REP. VEHÃC. AUTOM."
		replace sector = 8 if actividad_economica=="TRANSPORTE, ALMACENAMIENTO Y COMUNICACIONES"
		replace sector = 9 if actividad_economica=="HOTELES Y RESTAURANTES"
		replace sector = 10 if actividad_economica=="ACTIVIDADES INMOBILIARIAS, EMPRESARIALES Y DE ALQUILER"
		replace sector = 11 if actividad_economica=="INTERMEDIACIÃN FINANCIERA"
		replace sector = 12 if actividad_economica=="HOGARES PRIVADOS CON SERVICIO DOMÃSTICO"
		replace sector = 13 if actividad_economica=="ADMINISTRACIÃN PÃBLICA Y DEFENSA"
		replace sector = 14 if actividad_economica=="ENSEÃANZA"
		replace sector = 15 if actividad_economica=="SERVICIOS SOCIALES Y DE SALUD"
		replace sector = 16 if actividad_economica=="OTRAS ACTIV. SERV. COMUNITARIOS, SOCIALES Y PERSONALES"
		replace sector = 17 if actividad_economica=="ORGANIZACIONES Y ÃRGANOS EXTRATERRITORIALES"
		replace sector = 18 if actividad_economica=="NO DETERMINADO"
		label variable sector "Actividades económicas"
		label define sector 1 "Agropecuario" 2 "Pesca" 3 "Minería" 4 "Manufactura" 5 "Electricidad, Gas y Agua" 6 "Construcción" ///
							7 "Comercio" 8 "Transporte y comunicaciones" 9 "Hoteles y restaurantes" ///
							10 "Actividades inmobiliarias, empresariales y de alquiler" 11 "Intermediación financiera" ///
							12 "Servicio doméstico" 13 "Administración pública" 14 "Enseñanza" 15 "Servicios sociales y de salud" ///
							16 "Otras actividades" 17 "Organizaciones extraterritoriales" 18 "No determinado"
		label values sector sector

		** Sectores económicos agregado
		gen ss = .
		replace ss = 1 if inlist(sector,1)
		replace ss = 2 if inlist(sector,2)
		replace ss = 3 if inlist(sector,3)
		replace ss = 4 if inlist(sector,4)
		replace ss = 5 if inlist(sector,6)
		replace ss = 6 if inlist(sector,7)
		replace ss = 7 if ss==.
		label define ss 1 "Agropecuario" 2 "Pesca" 3 "Minería" 4 "Manufactura" 5 "Construcción" 6 "Comercio" 7 "Servicios"
		label values ss ss

		** Puestos de trabajo alrededor de la RMV
		gen rmv_5p = (remuneracion>=(1025/(1+0.05)) & remuneracion<=(1025*(1+0.05)))

	*** Tabulados

		** Trabajadores RMV
		table (tamaño), stat(mean rmv_5p)
		table (ss), stat(mean rmv_5p)

		** Información para el Ratio RMV (1025) / Salario promedio
		table (tamaño), stat(mean remuneracion)
		table (ss), stat(mean remuneracion)

	*** Gráfico

		** Información
		collapse (mean) rmv_5p, by(tamaño)
		replace rmv_5p = rmv_5p*100

		** Figura
		#delimit;
		graph hbar rmv_5p, 	over(tamaño, gap(*0.5) relabel(1 "Microempresa" 2 "Pequeña empresa" 3 "Mediana y Gran Empresa"))
							blabel(total, format(%2,1f))
							asyvar showyvars legend(off) bargap(20)
							nofill yscale(off) ytitle("")
							title("{bf:Trabajadores RMV}", span)
							bar(1, fintensity(100) fcolor(midblue%50) lcolor(black) lwidth(thin)) 
							bar(2, fintensity(100) fcolor(midblue%50) lcolor(black) lwidth(thin)) 
							bar(3, fintensity(100) fcolor(midblue%50) lcolor(black) lwidth(thin))
							text(30 20 "{bf:Total Formal}" "{bf:15,8%}", box size(medsmall) bcolor(none) 
																		  lcolor(black) lpattern(dash) bexpand bmargin(medsmall));
		#delimit cr
		graph export "${output}/graphs/G1.pdf", replace

**# Informalidad

	use "${data}/enaho01a-2022-500", clear

	*** Filtro del residente habitual
	keep if ((p204==1 & p205==2) | (p204==2 & p206==1)) & p501!=.

	*** Variables

		** Nivel educativo
		gen niv_educ = .
		replace niv_educ = 1 if p301a<=4 | p301a==12
		replace niv_educ = 2 if p301a==5 | p301a==6
		replace niv_educ = 3 if p301a==7 | p301a==8
		replace niv_educ = 4 if p301a==9 | p301a==10
		label define niv_educ 1 "Primaria" 2 "Secundaria" 3 "Superior no universitaria" 4 "Superior universitaria"
		label values niv_educ niv_educ

		** Empleo informal
		gen informal = 0 if ocu500==1
		replace informal = 1 if ocupinf==1
		label define informal 0 "Formal" 1 "Informal"
		label values informal informal
		
		** Ingreso laboral
		egen ingreso = rsum(i524a1 d529t i530a d536 i538a1 d540t i541a d543 d544t) if ocu500==1
		replace ingreso = ingreso/12

	*** Tabulados

		** Gráfico 3.1 - PEA ocupada con empleo formal e informal
		tab ocupinf if ocu500==1 [iw=fac500a]
		tab emplpsec if ocu500==1 [iw=fac500a]

		** Gráfico 3.5 - Tasa de empleo informal por características demográficas
		table (niv_educ) if ocu500==1 [iw=fac500a], stat(mean informal)

		** Gráfico 3.10 - Ingreso laboral real promedio mensual de la PEA ocupada con empleo formal e informal
		table (informal) if ocu500==1 & (ingreso>0 & !missing(ingreso)) [iw=fac500a], stat(mean ingreso) nformat(%10.0f)

**# Brecha de género

use "${data}/enaho01a-2022-500", clear

	*** Variables
	
		** Sexo
		gen sexo = 0
		replace sexo = 1 if p207==1
		lab var sexo "Sexo"
		lab define sexo 0 "Mujer" 1 "Hombre"
		lab values sexo sexo
	
		** Empleo informal
		gen informal = 0 if ocu500==1
		replace informal = 1 if ocupinf==1
		label define informal 0 "Formal" 1 "Informal"
		label values informal informal
		
		** Ingreso laboral
		egen ingreso = rsum(i524a1 d529t i530a d536 i538a1 d540t i541a d543 d544t) if ocu500==1
		replace ingreso = ingreso/12

		** Quintiles de ingreso
		xtile quintil = ingreso if ocu500==1 & (ingreso>0 & !missing(ingreso)) [aw=fac500a], nq(5)
		label define quintil 0 "Missing" 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4" 5 "Q5"
		label values quintil quintil

	*** Tabulados

		** Anexo 1.17 - PEA ocupada con empleo formal e informal
		table (sexo) if ocu500==1 & (ingreso>0 & !missing(ingreso)) [iw=fac500a], stat(mean ingreso) nformat(%10.0f)
		table (informal) (sexo) if ocu500==1 & (ingreso>0 & !missing(ingreso)) [iw=fac500a], stat(mean ingreso) nformat(%10.0f)		
		table (quintil) (sexo) if ocu500==1 & (ingreso>0 & !missing(ingreso)) [iw=fac500a], stat(mean ingreso) nformat(%10.0f)		

