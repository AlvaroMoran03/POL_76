' GeneraSmooth.prg
' por Diego Winkelried (Universidad del Pac�fico)
' Curso de Verano del BCRP - Econom�a Avanzada (Enero 2025)
' ----------------------------------------------------------------------------------------------------------------------
' Automatizaci�n del c�mputo de predicciones (en general, no �ptimas) recursivas y errores de predicci�n basadas en 
' m�todos de suavizamiento exponencial.
' ----------------------------------------------------------------------------------------------------------------------
' Este programa debe ser ejecutado despu�s en un Workfile que contenga la serie "y" y escalares "fin1", "fin2" y "Hmax"
' ----------------------------------------------------------------------------------------------------------------------
' Nota al instructor: Previamente, haga una demostraci�n del uso del "proc" SMOOTH
' ----------------------------------------------------------------------------------------------------------------------

' Par�metros del programa
' Etiqueta para los resultados (series y*, e* y grupos Gy*)
%tag = "s"

' Especificaci�n del m�todo de suavizamiento exponencial (S = "Simple", A = "Additive seasonality")
'%spec = "A"
'%freq = "12" ' Usar para la frecuencia en modelos con estacionalidad
%spec = "a" 
%freq = "12" 

' Limpieza de objetos posiblemente generados en ejecuciones anteriores de este programa
D y_* s_* y{%tag}* e{%tag}* Gy{%tag}*

' Muestras y horizontes
!T1 = fin1   ' Fin de la muestra de estimaci�n (fin1 es un escalar en el Workfile)
!T2 = fin2   ' Fin de la muestra total (fin2 es un escalar en el Workfile)
!maxH = Hmax ' M�ximo horizonte de proyecci�n (Hmax es un escalar en el Workfile)

' Estimaci�n del modelo de predicci�n
!T0 = @ROUND(!T2/4)
SMPL 1 !T0
y.SMOOTH({%spec}) ysm1_  {%freq}

!T0 = @ROUND(!T2/2)
SMPL 1 !T0
y.SMOOTH({%spec}) ysm2_ {%freq}

!T0 = @ROUND(3*!T2/4)
SMPL 1 !T0
y.SMOOTH({%spec}) ysm3_  {%freq}

SMPL @ALL
GROUP Gy{%tag} y ysm1_ ysm2_ ysm3_
FREEZE(Gy{%tag}_) Gy{%tag}.LINE 
SHOW Gy{%tag}_
D ysm*_ 

' Proyecciones recursivas: se generan las series y{%tag}1, y{%tag}2, ..., y{%tag}{!maxH} 
FOR !t = !T1+1 TO !T2-!maxH
	' Se proyecta en el horizonte t+1, t+2, ..., t+!maxH  
	SMPL 1 !t
	y.SMOOTH({%spec}) temp {%freq}

	' Guardamos los resultados para diferentes horizontes en diferentes series
	FOR !h = 1 TO !maxH 
		SMPL !t+!h !t+!h
		SERIES y{%tag}{!h} = temp
	NEXT !h

	D temp
NEXT !t

' Errores de predicci�n: se generan las series y{%tag}1, y{%tag}2, ..., y{%tag}{!maxH} 
SMPL @ALL
FOR !h = 1 TO !maxH 
	SERIES e{%tag}{!h} = y - y{%tag}{!h}
NEXT !h

' Creamos grupos para visualizar los resultados... y sus diferencias
SMPL @ALL

' Grupo 1a: Proyecciones para todos los horizontes, y la serie original
GROUP Gy{%tag}1a y{%tag}* y
Gy{%tag}1a.LINE
CLOSE Gy{%tag}1a

' Grupo 1b: Errores de predicci�n para todos los horizontes
GROUP Gy{%tag}1b e{%tag}* 0 
Gy{%tag}1b.LINE
CLOSE Gy{%tag}1b

' Grupo 2a: Proyecciones para horizontes seleccionados, y la serie original
GROUP Gy{%tag}2a y{%tag}1 y{%tag}{!maxH} y
Gy{%tag}2a.LINE
CLOSE Gy{%tag}2a

' Grupo 2b: Errores de predicci�n horizontes seleccionados
GROUP Gy{%tag}2b e{%tag}1 e{%tag}{!maxH} 0
Gy{%tag}2b.LINE
CLOSE Gy{%tag}2b

y.LINE
RUN EvaluacionOptima


