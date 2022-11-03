package Landwirtschaft
  package Hilfsmodelle
    model Motor_mech
      //Konstanten
      import Modelica.Constants.pi;
      //Parameter
      parameter Modelica.Units.SI.Inertia J = 0.5 * m * r ^ 2 "Massenträgheit des Motors";
      parameter Modelica.Units.SI.Mass m = 30 "Masse des Motors";
      parameter Modelica.Units.SI.Length r = 0.05 "Radius Motorwelle";
      parameter Real losses(final unit = "kg.m2/s") = 0.03 "Verlustkoeffizient Motor";
      parameter Modelica.Units.SI.Torque M_max = 25 "statisches Moment für Modus 1";
      //Variablen
      Modelica.Units.SI.AngularVelocity omega "Winkelgeschwindigkeit der Motorwelle";
      Modelica.Units.SI.Torque M_losses "Verlustmoment durch Reibung etc.";
      Modelica.Units.SI.Torque M_motor "Drehmoment des Motors";
      Modelica.Units.SI.Frequency n(displayUnit = "1/min") "Drehzahl des Motors";
      //Schnittstellen
      input Modelica.Blocks.Interfaces.RealInput vehicle_velocity(unit = "m/s") "aktuelle Geschwindigkeit des Traktors" annotation(
        Placement(visible = true, transformation(origin = {-80, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-59, -25}, extent = {{-15, -15}, {15, 15}}, rotation = 90)));
      input Modelica.Blocks.Interfaces.RealInput velocity_setPoint(unit = "m/s") "eingestellte Geschwindigkeit des Traktors" annotation(
        Placement(visible = true, transformation(origin = {-80, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-61, 27}, extent = {{-15, -15}, {15, 15}}, rotation = -90)));
      //Modelica-Blöcke
      Modelica.Blocks.Continuous.LimPID PowerController(yMax = M_max, yMin = 0, u_s = velocity_setPoint, u_m = vehicle_velocity, k = 600, Ti = 1, Td = 0.001, controllerType = Modelica.Blocks.Types.SimpleController.PID);
      Modelica.Blocks.Continuous.FirstOrder fo_unit(k = 1, T = 0.5);
      //Connectoren
      Landwirtschaft.Connectoren.Rotatoric abtrieb annotation(
        Placement(visible = true, transformation(origin = {82, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {101, -1}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    equation
      fo_unit.u = PowerController.y;
      M_motor = fo_unit.y;
      M_motor - M_losses - J * der(omega) = abtrieb.M;
      omega = der(abtrieb.alpha);
      n = omega / (2 * pi);
      M_losses = n ^ 2 * losses;
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(origin = {73, -1}, fillColor = {173, 173, 173}, fillPattern = FillPattern.HorizontalCylinder, extent = {{27, 9}, {-27, -9}}), Rectangle(origin = {-18, 1}, fillColor = {121, 121, 121}, fillPattern = FillPattern.Horizontal, extent = {{-64, 41}, {64, -41}}), Text(origin = {0, 50}, extent = {{-100, 7}, {100, -7}}, textString = "v_Setpoint"), Text(origin = {0, -47}, extent = {{-100, 7}, {100, -7}}, textString = "v_Measurement")}));
    end Motor_mech;

    model Getriebe
      //Parameter
      parameter Real i_a = 6.0 "Übersetzungsverhältnis Antrieb/Abtrieb";
      parameter Real i_z = 2.0 "Übersetzungsverhältnis Antrieb/Zapfwelle";
      parameter Real eta_a = 1.0 "Wirkungsgrad Antrieb/Abtrieb";
      parameter Real eta_z = 1.0 "Wirkungsgrad Antrieb/Zapfwelle";
      //Variablen
      Modelica.Units.SI.AngularVelocity omega_abtrieb;
      Modelica.Units.SI.AngularVelocity omega_antrieb;
      Modelica.Units.SI.AngularVelocity omega_zapfwelle;
      Modelica.Units.SI.AngularAcceleration omega_abtrieb_p;
      Modelica.Units.SI.AngularAcceleration omega_antrieb_p;
      Modelica.Units.SI.AngularAcceleration omega_zapfwelle_p;
      Modelica.Blocks.Interfaces.BooleanInput zapfwelle_Signal annotation(
        Placement(visible = true, transformation(origin = {-78, -74}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {38, 54}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      //Connectoren
      Landwirtschaft.Connectoren.Rotatoric antriebsseite annotation(
        Placement(visible = true, transformation(origin = {-88, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Connectoren.Rotatoric abtriebsseite annotation(
        Placement(visible = true, transformation(origin = {68, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {64, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b zapfwelle annotation(
        Placement(visible = true, transformation(origin = {66, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {62, 18}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
    equation
      if zapfwelle_Signal == true then
        omega_zapfwelle = omega_antrieb / i_z;
      else
        omega_zapfwelle = 0;
      end if;
      der(omega_abtrieb) = der(omega_antrieb) / i_a;
      abtriebsseite.M + antriebsseite.M * i_a * eta_a - zapfwelle.tau * i_z * eta_z = 0;
      omega_abtrieb = der(abtriebsseite.alpha);
      omega_antrieb = der(antriebsseite.alpha);
      omega_zapfwelle = der(zapfwelle.phi);
      omega_abtrieb_p = der(omega_abtrieb);
      omega_antrieb_p = der(omega_antrieb);
      omega_zapfwelle_p = der(omega_zapfwelle);
      annotation(
        Icon(graphics = {Rectangle(origin = {-24, 54}, fillColor = {113, 113, 113}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-54, 8}, {54, -8}}), Rectangle(origin = {10, -42}, fillColor = {113, 113, 113}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-54, 8}, {54, -8}}), Rectangle(origin = {-20, -35}, fillColor = {186, 186, 186}, fillPattern = FillPattern.Horizontal, extent = {{-10, 57}, {10, -57}}), Rectangle(origin = {-20, 55}, fillColor = {186, 186, 186}, fillPattern = FillPattern.Horizontal, extent = {{-14, 34}, {14, -34}}), Rectangle(origin = {11, 55}, fillColor = {186, 186, 186}, fillPattern = FillPattern.Horizontal, extent = {{-11, 13}, {11, -13}}), Rectangle(origin = {31, 18}, fillColor = {113, 113, 113}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-31, 4}, {31, -4}}), Line(origin = {59.71, 6.28}, points = {{-25, 0}, {5, 8.88178e-16}}, thickness = 1.75, arrow = {Arrow.Open, Arrow.Open}), Rectangle(origin = {11, 17}, fillColor = {186, 186, 186}, fillPattern = FillPattern.Horizontal, extent = {{-9, 25}, {9, -25}}), Text(origin = {-68, 77}, extent = {{-32, 13}, {32, -13}}, textString = "Antriebs-
    seite"), Text(origin = {45, -71}, extent = {{-53, 17}, {53, -17}}, textString = "Abtriebs-
    seite"), Text(origin = {64, 36}, extent = {{-34, 8}, {34, -8}}, textString = "Zapfwelle")}));
    end Getriebe;

    model Batterie
      //Teilmodelle
      Modelica.Electrical.Batteries.BatteryStacks.CellStack cellStack(cellData = cellData) annotation(
        Placement(visible = true, transformation(origin = {1, 0}, extent = {{-27, -20}, {27, 20}}, rotation = 0)));
      //Zelldaten für die Batteriezelle - Gibt die Ausgangsspannung in ABhängigkeit des Ladezustandes wieder.
      parameter Modelica.Electrical.Batteries.ParameterRecords.CellData cellData(Idis = 0.00001, OCV_SOC = [0.00, 0.735; 0.10, 0.882; 0.20, 0.946; 0.30, 0.956; 0.40, 0.963; 0.50, 0.965; 0.60, 0.967; 0.70, 0.971; 0.80, 0.976; 0.90, .0979; 1.00, 1.00], OCVmax = 27, OCVmin = 22, Qnom = 18000, R0 = cellData.Ri, Ri = cellData.OCVmax / 1200, SOCmax = 1, SOCmin = 0, T_ref = 293.15, alpha = 0, useLinearSOCDependency = true) annotation(
        Placement(visible = true, transformation(origin = {62, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //Connectoren
      Modelica.Electrical.Analog.Interfaces.PositivePin p annotation(
        Placement(visible = true, transformation(origin = {-89, -1}, extent = {{-27, -27}, {27, 27}}, rotation = 0), iconTransformation(origin = {70, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.NegativePin n annotation(
        Placement(visible = true, transformation(origin = {89, -1}, extent = {{-27, -27}, {27, 27}}, rotation = 0), iconTransformation(origin = {-72, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(p, cellStack.p) annotation(
        Line(points = {{-88, 0}, {-26, 0}}, color = {0, 0, 255}));
      connect(n, cellStack.n) annotation(
        Line(points = {{90, 0}, {28, 0}}, color = {0, 0, 255}));
      annotation(
        Icon(graphics = {Rectangle(origin = {-10, 0}, fillColor = {36, 198, 0}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-70, 40}, {70, -40}}), Rectangle(origin = {70, 0}, fillColor = {186, 186, 186}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{10, 20}, {-10, -20}}), Text(origin = {90, 30}, extent = {{-12, 30}, {12, -30}}, textString = "+"), Text(origin = {-91, 33}, extent = {{-11, 33}, {11, -33}}, textString = "-")}));
    end Batterie;

    model Rad
      //Parameter
      parameter Modelica.Units.SI.Mass m = 20 "Masse";
      parameter Modelica.Units.SI.Inertia J = 0.5 * m * r ^ 2 "Massenträgheit";
      parameter Modelica.Units.SI.Length r = 0.6 "Radius Radmitte - Rollfläche";
      //Variablen
      Modelica.Units.SI.AngularVelocity omega "Winkelgeschwindigkeit Rad";
      //Connectoren
      Landwirtschaft.Connectoren.Translatoric translatoric annotation(
        Placement(visible = true, transformation(origin = {0, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //Schnittstellen
      Modelica.Blocks.Interfaces.RealOutput vehicle_vel annotation(
        Placement(visible = true, transformation(origin = {-16, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-102, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Landwirtschaft.Connectoren.Rotatoric rotatoric annotation(
        Placement(visible = true, iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      vehicle_vel = der(translatoric.s);
      J * der(omega) = rotatoric.M - translatoric.F * r;
      omega = der(rotatoric.alpha);
      rotatoric.alpha = translatoric.s / r;
      annotation(
        Icon(graphics = {Text(origin = {-43, 60}, extent = {{-57, 12}, {57, -12}}, textString = "Radlager"), Text(origin = {-45, 1}, extent = {{-53, 11}, {53, -11}}, textString = "Radwelle"), Text(origin = {-27, -56}, extent = {{-59, 10}, {59, -10}}, textString = "v_Measurment"), Rectangle(extent = {{-100, 100}, {100, -100}}), Ellipse(origin = {70, 70}, fillPattern = FillPattern.Solid, extent = {{-30, 30}, {30, -30}}), Ellipse(origin = {70, 70}, fillColor = {211, 211, 211}, fillPattern = FillPattern.Solid, extent = {{-22, 22}, {22, -22}}), Ellipse(origin = {70, 86}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {86, 70}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {70, 54}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {54, 70}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {82, 58}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {60, 58}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {60, 82}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {82, 82}, fillPattern = FillPattern.Solid, extent = {{-2, 2}, {2, -2}}), Ellipse(origin = {70, 70}, fillPattern = FillPattern.Solid, extent = {{-4, 4}, {4, -4}})}));
    end Rad;

    model Chassis
      //Parameter
      parameter Modelica.Units.SI.Mass m = 400 "Gewicht Traktor";
      parameter Real cW = 0.95 "Luftwiederstandsbeiwert";
      parameter Modelica.Units.SI.Area A = 5 "Fläche Traktor vorne";
      constant Modelica.Units.SI.Density rho = 1.204 "Dichte von Luft bei 20°C und 1,013 bar";
      //Variablen
      Modelica.Units.SI.Acceleration a "Translatorische Beschleunigung in horizontale";
      Modelica.Units.SI.Velocity v "Translatorische Geschwindigkeit in horizontale";
    public
      Modelica.Units.SI.Force F_air "Luftwiderstand bei der Bewegung";
      //Connectoren
      Landwirtschaft.Connectoren.Translatoric translatoric annotation(
        Placement(visible = true, transformation(origin = {40, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {52, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Connectoren.Translatoric anhaengerkupplung annotation(
        Placement(visible = true, transformation(origin = {84, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      F_air = 0.5 * cW * A * rho * v ^ 2;
      m * a = translatoric.F + anhaengerkupplung.F - F_air;
      der(v) = a;
      der(translatoric.s) = v;
      translatoric.s = anhaengerkupplung.s;
      annotation(
        Icon(coordinateSystem(extent = {{-100, -50}, {80, 50}}), graphics = {Rectangle(origin = {-10, 0}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, extent = {{-90, 50}, {90, -50}}), Text(origin = {48, 30}, rotation = 180, extent = {{-30, 18}, {30, -18}}, textString = "Anhänger-
kupplung"), Text(origin = {49, -34}, extent = {{-27, 8}, {27, -8}}, textString = "Radlager"), Text(origin = {-48, 24}, extent = {{-50, 24}, {50, -24}}, textString = "Masse: %m
Fläche Front: %A
cW-Wert: %cW")}),
        Diagram(coordinateSystem(extent = {{-100, -50}, {80, 50}})));
    end Chassis;

    model Kegelradgetriebe
      // Konstanten
      //Parameter
      parameter Real i_ka = 2.0 "Übersetzungsverhältnis Antrieb/Abtrieb";
      parameter Real eta_k = 1.0 "Wirkungsgrad Antrieb/Abtrieb";
      //Variablen
      Modelica.Units.SI.AngularVelocity omegaK_abtrieb;
      Modelica.Units.SI.AngularVelocity omegaK_antrieb;
      Modelica.Units.SI.AngularAcceleration omegaK_abtrieb_a;
      Modelica.Units.SI.AngularAcceleration omegaK_antrieb_a;
      // Visual
      //Gleichungen
      Modelica.Mechanics.Rotational.Interfaces.Flange_a K_Antriebsseite annotation(
        Placement(visible = true, transformation(origin = {-52, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-82, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b K_Abtriebsseite annotation(
        Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {18, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      der(omegaK_abtrieb) = der(omegaK_antrieb) / i_ka;
      omegaK_abtrieb_a = der(omegaK_abtrieb);
      omegaK_antrieb_a = der(omegaK_antrieb);
      omegaK_abtrieb = der(K_Abtriebsseite.phi);
      omegaK_antrieb = der(K_Antriebsseite.phi);
      K_Abtriebsseite.tau - K_Antriebsseite.tau * i_ka * eta_k = 0;
//Visual
      annotation(
        Icon(graphics = {Polygon(origin = {-10, -40}, fillColor = {148, 148, 148}, fillPattern = FillPattern.Forward, lineThickness = 1, points = {{-10, 40}, {-10, -40}, {10, -20}, {10, 20}, {-10, 40}, {-10, 40}}), Rectangle(origin = {-50, -39}, fillColor = {98, 98, 98}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{30, -11}, {-30, 11}}), Rectangle(origin = {19, 30}, fillColor = {98, 98, 98}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-11, 30}, {11, -30}}), Polygon(origin = {20, -10}, fillColor = {147, 147, 147}, fillPattern = FillPattern.Forward, lineThickness = 1, points = {{-20, -10}, {20, -10}, {40, 10}, {-40, 10}, {-40, 10}, {-20, -10}}), Line(origin = {20, -40}, points = {{-40, -40}, {40, 40}, {40, 40}}, pattern = LinePattern.Dash, thickness = 0.5), Line(origin = {0, -20}, points = {{-20, 20}, {20, -20}, {20, -20}}, pattern = LinePattern.Dash, thickness = 0.5), Text(origin = {-53, -21}, extent = {{-19, 5}, {19, -5}}, textString = "Antrieb"), Text(origin = {47, 7}, extent = {{-19, 5}, {19, -5}}, textString = "Abtrieb")}));
    end Kegelradgetriebe;

    model Drehscheibe
      // Parameter
      parameter Modelica.Units.SI.Mass m = 50 "Gewicht der Drehscheibe";
      parameter Modelica.Units.SI.Inertia J_D = 20 "Trägheitsmoment der Drehscheibe";
      // Variablen
      Modelica.Units.SI.AngularVelocity omegaD_abtrieb;
      Modelica.Units.SI.AngularAcceleration omegaD_abtrieb_a;
  Modelica.Mechanics.Rotational.Interfaces.Flange_a D_Abtriebsseite annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      omegaD_abtrieb = der(D_Abtriebsseite.phi);
      omegaD_abtrieb_a = der(omegaD_abtrieb);
      omegaD_abtrieb_a * J_D + D_Abtriebsseite.tau = 0;
      annotation(
        Icon(graphics = {Ellipse(origin = {0, -1}, fillColor = {104, 104, 104}, fillPattern = FillPattern.CrossDiag, extent = {{-54, 53}, {54, -53}})}));
    end Drehscheibe;

    model Zahnrad
      // Parameter
      parameter Real eta_z = 0.7 "Wirkungsgrad Antrieb/Abtrieb";
      parameter Real i_z = -1.0;
      // Variablen
      Modelica.Units.SI.AngularVelocity omegaZ_abtrieb;
      Modelica.Units.SI.AngularVelocity omegaZ_antrieb;
      Modelica.Units.SI.AngularAcceleration omegaZ_abtrieb_a;
      Modelica.Units.SI.AngularAcceleration omegaZ_antrieb_a;
        Modelica.Mechanics.Rotational.Interfaces.Flange_a Z_Antriebsseite annotation(
        Placement(visible = true, transformation(origin = {-54, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-20, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b Z_Abtriebsseite annotation(
        Placement(visible = true, transformation(origin = {46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {20, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Gleichungen
    equation
      der(omegaZ_abtrieb) = der(omegaZ_antrieb) / i_z;
      omegaZ_antrieb_a = der(omegaZ_antrieb);
      omegaZ_abtrieb_a = der(omegaZ_abtrieb);
      omegaZ_antrieb = der(Z_Antriebsseite.phi);
      omegaZ_abtrieb = der(Z_Abtriebsseite.phi);
      Z_Antriebsseite.tau + Z_Abtriebsseite.tau * eta_z = 0;
// Darstellung
      annotation(
        Icon(graphics = {Rectangle(fillColor = {121, 121, 121}, fillPattern = FillPattern.Horizontal, lineThickness = 1, extent = {{-20, 80}, {20, -80}}), Text(origin = {-36, 13}, extent = {{-14, 15}, {14, -15}}, textString = "Antrieb"), Text(origin = {36, 13}, extent = {{-14, 15}, {14, -15}}, textString = "Abtrieb")}));
    end Zahnrad;

    model Tank
    // Konstanten
    constant Modelica.Units.SI.Acceleration g = Modelica.Constants.g_n;
    constant Real roh = 1000 "Dichte von Wasser";
    // Parameter
    parameter Modelica.Units.SI.Pressure p_e = 1.013*10^5 "Umgebungsdruck in bar";
    parameter Modelica.Units.SI.Area A = 0.1 "Fläche der Zylinders";
    parameter Modelica.Units.SI.Area At = 0.01 "Fläche der Öffnung"; 
// Variablen
    Modelica.Units.SI.Distance h "Höhendifferenz";
    Modelica.Units.SI.Pressure p_a "Umgebungsdruck";
    Modelica.Units.SI.Velocity v_t "Geschwindigkeit Tank";
    Modelica.Units.SI.Acceleration a_t "Beschleunigung Tank";
    Real V (start = 1) "Tankvolumen";
    Real c;
    Real V_t;
    Real m;
    
    
    Connectoren.Translatoric translatoric_tank annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    
    // Gleichungen
equation
// Druckausgleich
      p_e = p_a;
// Austrittsgeschwindigkeit
      c = sqrt(2 * g * h);
// Füllvolumen
      V = A * h;
// Volumenstrom
      der(V) = -V_t;
// Volumenstrom
      V_t = c * At;
    // Masse Tank
      translatoric_tank.F = V * roh * g;
    // Weg Tank
      a_t = der(v_t);
      v_t = der(translatoric_tank.s)
    
    annotation(
        Icon(graphics = {Rectangle(origin = {0, 20}, fillPattern = FillPattern.Solid, extent = {{-40, -60}, {40, 60}})}));
    
    end Tank;
    
    model Tank_Test
    // Konstanten
    constant Modelica.Units.SI.Acceleration g = Modelica.Constants.g_n;
    constant Real roh = 1000 "Dichte von Wasser";
    // Parameter
    parameter Modelica.Units.SI.Pressure p_e = 1.013*10^5 "Umgebungsdruck in bar";
    parameter Modelica.Units.SI.Area A = 0.1 "Fläche der Zylinders";
    parameter Modelica.Units.SI.Area At = 0.01 "Fläche der Öffnung"; 
// Variablen
    Modelica.Units.SI.Distance h "Höhendifferenz";
    Modelica.Units.SI.Pressure p_a "Umgebungsdruck";
    Real V (start = 1) "Tankvolumen";
    Real c;
    Real V_t;
    // Gleichungen

equation
// Druckausgleich
      p_e = p_a;
// Austrittsgeschwindigkeit
      c = sqrt(2 * g * h);
// Füllvolumen
      V = A * h;
// Volumenstrom
      der(V) = -V_t;
// Volumenstrom
      V_t = c * At;
    annotation(
        Icon(graphics = {Rectangle(origin = {0, 20}, fillPattern = FillPattern.Solid, extent = {{-40, -60}, {40, 60}})}));end Tank_Test;
  end Hilfsmodelle;

  package Connectoren
    connector Translatoric
      flow Modelica.Units.SI.Force F;
      Modelica.Units.SI.Distance s;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {10, 149, 0}, fillPattern = FillPattern.Solid, extent = {{-60, 60}, {60, -60}})}));
    end Translatoric;

    connector Rotatoric
      flow Modelica.Units.SI.Torque M;
      Modelica.Units.SI.Angle alpha;
      annotation(
        Icon(graphics = {Ellipse(fillColor = {197, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-60, 60}, {60, -60}})}));
    end Rotatoric;

    connector Fluid
    
      Modelica.Units.SI.Pressure w_t;           // Strömujngsgeschwindigkeit
      flow Real V_t;                      // Volumenstrom
    
      annotation(
        Icon(graphics = {Ellipse(origin = {0, 3}, fillColor = {48, 133, 190}, fillPattern = FillPattern.Solid, extent = {{-64, 65}, {64, -65}})}));
    end Fluid;
  end Connectoren;

  package Arbeitsgeraete
  model Masse
      // Konstanten
      // Parameter
      parameter Modelica.Units.SI.Mass m = 5000 "Masse in kg";
      // Variable
      Modelica.Units.SI.Acceleration a_m;
      Modelica.Units.SI.Velocity v_m;
      Landwirtschaft.Connectoren.Translatoric translatoric annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {2, 0}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
      // Gleichungen
    equation
      der(translatoric.s) = v_m;
      der(v_m) = a_m;
      translatoric.F = m * a_m;
// Beschreibung
      annotation(
        Icon(graphics = {Polygon(origin = {5, -7}, fillColor = {108, 108, 108}, fillPattern = FillPattern.Solid, points = {{-71, 41}, {-33, 69}, {57, 55}, {71, 9}, {59, -65}, {-23, -43}, {-81, -45}, {-61, 3}, {-61, 3}, {-71, 41}})}));
    end Masse;

    model Duengerverteiler
    Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {-31, -41}, extent = {{-45, -45}, {45, 45}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Zahnrad zahnrad annotation(
        Placement(visible = true, transformation(origin = {37, 1}, extent = {{-37, -37}, {37, 37}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe1 annotation(
        Placement(visible = true, transformation(origin = {-16, 66}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe2 annotation(
        Placement(visible = true, transformation(origin = {16, 66}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
    equation
      connect(kegelradgetriebe.K_Abtriebsseite, drehscheibe1.D_Abtriebsseite) annotation(
        Line(points = {{-22, -14}, {-60, -14}, {-60, 66}, {-16, 66}}));
  connect(kegelradgetriebe.K_Abtriebsseite, zahnrad.Z_Antriebsseite) annotation(
        Line(points = {{-22, -14}, {-6, -14}, {-6, 0}, {30, 0}}));
  connect(zahnrad.Z_Abtriebsseite, drehscheibe2.D_Abtriebsseite) annotation(
        Line(points = {{44, 0}, {80, 0}, {80, 66}, {16, 66}}));
    annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "D", fontSize = 150)}));end Duengerverteiler;

    model Anhaenger
    
        //Parameter
        parameter Modelica.Units.SI.Mass a_masse = 2895 "Anhänger Masse in kg";
        parameter Modelica.Units.SI.Mass gr_masse = 5500 "Granulat Masse in kg";
        parameter Modelica.Units.SI.Mass Ent_masse = 10 "Entlassungsmasse der Granulaten/sek";
        parameter Real cW = 0.95 "Luftwiederstandsbeiwert";
        parameter Modelica.Units.SI.Area A_F = 15 "Fläche Anhänger";
        parameter Modelica.Units.SI.Area T_F = 10 "Fläche Traktor";
        constant Modelica.Units.SI.Density dichte = 1.204 "Dichte von Luft bei 20°C und 1,013 bar";
        
        //Variable
        Modelica.Units.SI.Acceleration a "Translatorische Beschleunigung in horizontale";
        Modelica.Units.SI.Velocity v "Translatorische Geschwindigkeit in horizontale";
        Modelica.Units.SI.Force F_Luft "Luftwiderstand bei der Bewegung";
        
        //Klasse
        Arbeitsgeraete.AnhaengerZustand zustand(anhaenger_masse = a_masse,granulat_masse = gr_masse,entlassungs_masse = Ent_masse,     
        istbefuellt = true);
        
        //Connectoren
        Landwirtschaft.Connectoren.Translatoric anhaenger_kupplung annotation(
          Placement(visible = true, transformation(extent = {{0, 0}, {0, 0}}, rotation = 0), iconTransformation(origin = {-91, 11}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
        Landwirtschaft.Connectoren.Translatoric translatoric_radlager annotation(
          Placement(visible = true, transformation(origin = {46, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {46, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
     
        // Gleichungen
        equation
        F_Luft = 0.5 * cW * A_F * dichte * v ^ 2;
        zustand.gesamte_anhaenger_masse * a = translatoric_radlager.F + anhaenger_kupplung.F - F_Luft;
        der(v) = a;
        der(anhaenger_kupplung.s) = v;
        translatoric_radlager.s = anhaenger_kupplung.s;
        //Beschreibung
        annotation(
            Icon(graphics = {Rectangle(origin = {-67, 10}, extent = {{-25, 2}, {25, -2}}), Ellipse(origin = {46, -13}, pattern = LinePattern.Dash, extent = {{-14, 15}, {14, -15}}), Ellipse(origin = {48, 33}, pattern = LinePattern.Dash, extent = {{-14, 15}, {14, -15}}), Ellipse(origin = {-20, -11}, pattern = LinePattern.Dash, extent = {{-14, 15}, {14, -15}}), Polygon(origin = {14, 7}, fillColor = {98, 98, 98}, pattern = LinePattern.Dash, fillPattern = FillPattern.Horizontal, points = {{-50, 19}, {48, 19}, {70, 39}, {70, -41}, {48, -21}, {-50, -21}, {-70, -41}, {-70, 41}, {-50, 19}}), Ellipse(origin = {-20, 33}, pattern = LinePattern.Dash, extent = {{-14, 15}, {14, -15}}), Polygon(origin = {22, 12}, fillColor = {24, 24, 24}, fillPattern = FillPattern.Solid, points = {{-34, 30}, {74, 30}, {38, -30}, {-64, -30}, {-66, 30}, {-66, 30}, {-34, 30}}), Text(origin = {14, 15}, lineColor = {255, 255, 255}, extent = {{-32, 7}, {32, -7}}, textString = "Behälter"), Text(origin = {83, -8}, lineColor = {255, 0, 0}, extent = {{-23, 6}, {23, -6}}, textString = "chassis"), Text(origin = {14, 15}, lineColor = {255, 255, 255}, extent = {{-32, 7}, {32, -7}}, textString = "Behälter"), Text(origin = {14, 15}, lineColor = {255, 255, 255}, extent = {{-32, 7}, {32, -7}}, textString = "Behälter"), Text(origin = {50, -31}, extent = {{-12, -5}, {12, 5}}, textString = "Radlager"), Text(origin = {-90, -1}, extent = {{-32, 7}, {32, -7}}, textString = "Anhängerkopplung")}));
    end Anhaenger;

    model Anhaenger_Rad
    
      //Konstanten
      constant Modelica.Units.SI.Acceleration g =Modelica.Constants.g_n "Erdbeschleunigung";
      constant Real reibungskoeff = 0.45 "Reibungskoeffizient sowohl auf dem Erdweg als auch auf dem Ackerboden";
      
      //Parameter
      parameter Modelica.Units.SI.Mass m = 284.36 "Rad Masse";
      parameter Modelica.Units.SI.Inertia J = 0.5 * m * r ^ 2 "Massenträgheit";
      parameter Modelica.Units.SI.Length r = 0.8985 "Radius Radmitte - Rollfläche";
      
      //Variablen
      Modelica.Units.SI.AngularVelocity omega "Winkelgeschwindigkeit Rad";
      Modelica.Units.SI.Force F_roll "Rollreibungskraft";
      
      //Connectoren
      Landwirtschaft.Connectoren.Translatoric radlager annotation(
        Placement(visible = true, transformation(origin = {0, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput vehicle_vel annotation(
        Placement(visible = true, transformation(origin = {-16, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-74, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Landwirtschaft.Connectoren.Rotatoric radwelle annotation(
        Placement(visible = true, transformation(extent = {{0, 0}, {0, 0}}, rotation = 0), iconTransformation(origin = {-80, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      
      equation
      vehicle_vel = der(radlager.s);
      //F_roll = reibungskoeff * m * g;
      if vehicle_vel > 0 then
        F_roll = reibungskoeff * m * g;
      else
        F_roll = 0;
      end if;
      J * der(omega) = radwelle.M - radlager.F * r - F_roll*r;
      omega = der(radwelle.alpha);
      radwelle.alpha = radlager.s / r;
     
      annotation(
        Icon(graphics = {Ellipse(origin = {4, 16}, fillPattern = FillPattern.Solid, extent = {{-38, 38}, {38, -38}}), Ellipse(origin = {5, 16}, fillColor = {134, 134, 134}, fillPattern = FillPattern.Solid, extent = {{29, -30}, {-29, 30}}), Polygon(origin = {-13, -12}, fillColor = {139, 139, 139}, fillPattern = FillPattern.Solid, points = {{11, 30}, {-19, -24}, {-9, -32}, {23, 24}, {11, 30}}), Ellipse(origin = {-26, -35}, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid, extent = {{-6, 17}, {6, -17}}), Ellipse(origin = {4, 15}, fillColor = {115, 115, 115}, fillPattern = FillPattern.Solid, extent = {{-6, 17}, {6, -17}}), Line(origin = {5, 18}, points = {{-29, 0}, {29, 0}}, color = {255, 255, 255}), Line(origin = {5, 16}, points = {{1, 30}, {-1, -30}}, color = {255, 255, 255}), Line(origin = {5, 17}, points = {{-19, 23}, {19, -23}}, color = {255, 255, 255}), Line(origin = {4, 17}, points = {{22, 19}, {-22, -19}}, color = {255, 255, 255}), Text(origin = {-54, -9}, extent = {{-24, 5}, {24, -5}}, textString = "Radwelle"), Text(origin = {-55, 10}, extent = {{-15, 6}, {15, -6}}, textString = "Radlager")}));
    end Anhaenger_Rad;

    class AnhaengerZustand
    
      //Parameter
      parameter Modelica.Units.SI.Mass anhaenger_masse = 2895 "Anhänger Masse in kg";
      parameter Modelica.Units.SI.Mass granulat_masse = 5500 "Granulat Masse in kg";
      parameter Modelica.Units.SI.Mass entlassungs_masse = 10 "Entlassungsmasse der Granulaten/sek";
      parameter Boolean istbefuellt = true "gibt an ob der Anhänger mit Granulaten befüllt ist";
    
      //Variable
      Modelica.Units.SI.Mass gesamte_anhaenger_masse;
      Real b "Konstante b der linearen Funktion" ;
      Real k "Steigung der linearen Funktion"; // m = k*t + b
      Real gesamte_zeit "notwendige Zeit bis alle Granulate entlassen werden";
      Real t(start=0) "time";
    
      //Gleichungen
      equation
      if istbefuellt then
       gesamte_zeit = granulat_masse / entlassungs_masse;
       b = anhaenger_masse + granulat_masse;
       k = (anhaenger_masse - b) / gesamte_zeit;
       der(t) = if t>=gesamte_zeit then 0 else 1;
       gesamte_anhaenger_masse = k*t+b;
      else
        gesamte_anhaenger_masse = anhaenger_masse;
      end if;
     
      end AnhaengerZustand;
  end Arbeitsgeraete;

  package Zugmaschinen
    model Traktor
      //Parametrierung
      //Motor
      parameter Modelica.Units.SI.Mass m_Motor = 30.0 "Masse der Motorwelle";
      parameter Modelica.Units.SI.Length r_Motor = 0.05 "Radius der Motorwelle";
      parameter Real loss_Motor(final unit = "kg.m2/s") = 0.03 "Motorverluste in Abhängigkeit von Drehzahl";
      parameter Modelica.Units.SI.Torque M_Motor_max = 380 "Nennmoment des Motors";
      //Getriebe
      parameter Real i_Getriebe(final unit = "1") = 6 "Übersetzungsverhältnis des Getriebes";
      parameter Real eta_Getriebe(final unit = "1") = 0.98 "Wirkungsgrad des Getriebes";
      parameter Real i_Zapfwelle(final unit = "1") = 2 "Übersetzungsverhältnis des Getriebes";
      parameter Real eta_Zapfwelle(final unit = "1") = 0.95 "Wirkungsgrad des Getriebes";
      //Rad
      parameter Modelica.Units.SI.Length r_Rad = 0.9 "Radius des angetriebenen Rades";
      parameter Modelica.Units.SI.Mass m_Rad = 400 "Masse des Rades";
      //Chassis
      parameter Modelica.Units.SI.Mass m_Chassis = 1500 "Masse des Chassis";
      parameter Real cW_Chassis = 0.95 "Luftwiederstandsbeiwert";
      parameter Modelica.Units.SI.Area A_Chassis = 5 "Fläche Traktor vorne";
      //Teilmodelle - Eigen
      Landwirtschaft.Hilfsmodelle.Chassis chassis(m = m_Chassis, cW = cW_Chassis, A = A_Chassis) annotation(
        Placement(visible = true, transformation(origin = {135.841, 223.584}, extent = {{-108.832, -54.4158}, {87.0652, 54.4158}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Motor_mech motor_mech(M_max = M_Motor_max, losses = loss_Motor, r = r_Motor, m = m_Motor) annotation(
        Placement(visible = true, transformation(origin = {-163, 99}, extent = {{-47, -47}, {47, 47}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Getriebe getriebe(i_a = i_Getriebe, i_z = i_Zapfwelle, eta_a = eta_Getriebe, eta_z = eta_Zapfwelle) annotation(
        Placement(visible = true, transformation(origin = {-38, 64}, extent = {{-62, -62}, {62, 62}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Rad rad(m = m_Rad, r = r_Rad) annotation(
        Placement(visible = true, transformation(origin = {178, 16}, extent = {{-72, -72}, {72, 72}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Batterie batterie annotation(
        Placement(visible = true, transformation(origin = {-157, -199}, extent = {{-43, -43}, {43, 43}}, rotation = 90)));
      //Teilmodelle - Modelica
      Modelica.Electrical.Analog.Ideal.IdealClosingSwitch switch annotation(
        Placement(visible = true, transformation(origin = {-2, -146}, extent = {{-34, -34}, {34, 34}}, rotation = 0)));
      Modelica.Electrical.Analog.Basic.Ground ground annotation(
        Placement(visible = true, transformation(origin = {13, -283}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
      //Connectoren - Eigen
      Landwirtschaft.Connectoren.Translatoric anhaengerkupplung annotation(
        Placement(visible = true, transformation(origin = {276, 224}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {249, 69}, extent = {{-25, -25}, {25, 25}}, rotation = 0)));
      //Connectoren - Modelica
      Modelica.Electrical.Analog.Interfaces.PositivePin p annotation(
        Placement(visible = true, transformation(origin = {280, -258}, extent = {{-12, -12}, {12, 12}}, rotation = 0), iconTransformation(origin = {249, -15}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      Modelica.Electrical.Analog.Interfaces.NegativePin n annotation(
        Placement(visible = true, transformation(origin = {281, -145}, extent = {{-15, -15}, {15, 15}}, rotation = 0), iconTransformation(origin = {249, 9}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a zapfwelle annotation(
        Placement(visible = true, transformation(origin = {276, 102}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {250, 38}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
      //Eingangs-/ Steuersignale
      Modelica.Blocks.Interfaces.RealInput vel_Setpoint annotation(
        Placement(visible = true, transformation(origin = {-300, 140}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-280, 202}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.BooleanInput Zapfwellensignal annotation(
        Placement(visible = true, transformation(origin = {-300, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-280, 120}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.BooleanInput Batterieschalter annotation(
        Placement(visible = true, transformation(origin = {-300, -130}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-280, 76}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.BooleanOutput Steuersignal annotation(
        Placement(visible = true, transformation(origin = {293, -279}, extent = {{-13, -13}, {13, 13}}, rotation = 0), iconTransformation(origin = {248, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.BooleanInput Steuersingal_Input annotation(
        Placement(visible = true, transformation(origin = {-300, -280}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-280, 32}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    equation
      connect(rad.vehicle_vel, motor_mech.vehicle_velocity) annotation(
        Line(points = {{105, -27}, {-189.12, -27}, {-189.12, 88.8}}, color = {0, 0, 127}));
      connect(rad.rotatoric, getriebe.abtriebsseite) annotation(
        Line(points = {{106, 16}, {2, 16}, {2, 38}}));
      connect(getriebe.antriebsseite, motor_mech.abtrieb) annotation(
        Line(points = {{-86.36, 97.48}, {-116.36, 97.48}}));
      connect(rad.translatoric, chassis.translatoric) annotation(
        Line(points = {{106, 59}, {183, 59}, {183, 169}}));
      connect(switch.p, batterie.p) annotation(
        Line(points = {{-36, -146}, {-156, -146}, {-156, -169}}, color = {0, 0, 255}));
      connect(anhaengerkupplung, chassis.anhaengerkupplung) annotation(
        Line(points = {{276, 224}, {213, 224}}));
      connect(switch.n, n) annotation(
        Line(points = {{32, -146}, {281, -146}, {281, -145}}, color = {0, 0, 255}));
      connect(vel_Setpoint, motor_mech.velocity_setPoint) annotation(
        Line(points = {{-300, 140}, {-192, 140}, {-192, 112}}, color = {0, 0, 127}));
      connect(Zapfwellensignal, getriebe.zapfwelle_Signal) annotation(
        Line(points = {{-300, -40}, {40, -40}, {40, 98}, {-14, 98}}, color = {255, 0, 255}));
      connect(Batterieschalter, switch.control) annotation(
        Line(points = {{-300, -130}, {-58, -130}, {-58, -84}, {-2, -84}, {-2, -106}}, color = {255, 0, 255}));
      connect(getriebe.zapfwelle, zapfwelle) annotation(
        Line(points = {{0, 76}, {72, 76}, {72, 102}, {276, 102}}));
      connect(Steuersingal_Input, Steuersignal) annotation(
        Line(points = {{-300, -280}, {294, -280}, {294, -278}}, color = {255, 0, 255}));
      connect(ground.p, p) annotation(
        Line(points = {{13, -258}, {280, -258}}, color = {0, 0, 255}));
      connect(ground.p, batterie.n) annotation(
        Line(points = {{13, -258}, {-157, -258}, {-157, -230}}, color = {0, 0, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-300, -300}, {300, 300}}), graphics = {Rectangle(origin = {12, 114}, fillColor = {247, 238, 139}, fillPattern = FillPattern.Solid, extent = {{-252, 184}, {252, -184}}), Text(origin = {-148, 229}, extent = {{-88, 19}, {88, -19}}, textString = "Mechanical"), Rectangle(origin = {12, -184}, fillColor = {170, 203, 255}, fillPattern = FillPattern.Solid, extent = {{-252, 114}, {252, -114}}), Text(origin = {-151, -104}, extent = {{-81, 18}, {81, -18}}, textString = "Electrical")}),
        Icon(coordinateSystem(extent = {{-300, -300}, {300, 300}}), graphics = {Bitmap(origin = {13, -67}, extent = {{-297, -195}, {297, 195}}, fileName = "modelica://Landwirtschaft/Traktor.png"), Text(origin = {-190, 208}, extent = {{-70, 22}, {70, -22}}, textString = "vel_Setpoint"), Text(origin = {-202, 120}, extent = {{-70, 16}, {70, -16}}, textString = "PTO-Shaft"), Text(origin = {-176, 76}, extent = {{-102, 16}, {102, -16}}, textString = "Battery-Switch"), Text(origin = {-180, 34}, extent = {{-102, 16}, {102, -16}}, textString = "Control-Signal")}),
        Documentation(info = "<html><head></head><body><div>Das Traktormodell besteht aus einer Menge von Submodellen und Komponenten:</div><div><ul><li>7 Teilmodelle:</li><ul><li>Motor</li><li>Getriebe</li><li>Rad</li><li>Chassis</li></ul><li>4 Connectoren zur Anbindung von Arbeitsgeräten:</li><ul><li>mechanisch translatorische Schnittstelle (repräsentiert die Anhängerkupplung)</li><li>mechanisch rotatorische Schnittstelle (repräsentiert die Zapfwelle)</li><li>2 elektrische Schnittstellen (repräsentieren den + und – Pol)</li></ul><li>4 Eingangssignale zur Steuerung des Traktorverhaltens:</li><ul><li>Real-Input zur Steuerung der Zielgeschwindigkeit</li><li>Boolscher Wert zur Steuerung des eingelegten Ganges</li><li>Boolscher Wert zur Steuerung der Zapfwelle</li><li>Boolscher Wert zur Steuerung des Stromkreises der Batterie</li></ul></ul></div><div>Im Folgenden soll kurz auf die Parameter des Traktors eingegangen werden. Neben den 4 Eingangsparametern, besitzt der Traktor 13 Parameter, welche direkt in dessen Parametrierungsfenster eingestellt werden können und entsprechend an die Submodelle weitergeleitet werden:</div><div><br></div><div><b>Beschreibung<span class=\"Apple-tab-span\" style=\"white-space:pre\">						</span>Variablenname<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Standardwert&nbsp;<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>Einheit</b></div><div>Masse der Motorwelle<span class=\"Apple-tab-span\" style=\"white-space:pre\">					</span>m_Motor<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span><span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>30.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>kg</div><div>Radius der Motorwelle<span class=\"Apple-tab-span\" style=\"white-space:pre\">					</span>r_Motor<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>0.05<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>m</div><div>Motorverluste<span class=\"Apple-tab-span\" style=\"white-space:pre\">						</span>loss_Motor<span class=\"Apple-tab-span\" style=\"white-space: pre;\">		</span><span class=\"Apple-tab-span\" style=\"white-space: pre;\">	</span>0.03<span class=\"Apple-tab-span\" style=\"white-space: pre;\">			</span>Nm/s</div><div>maximales Moment<span class=\"Apple-tab-span\" style=\"white-space:pre\">					</span>M_Motor_max<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>380.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>Nm</div><div>Übersetzung Getriebe AntriebAbtrieb<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>i_Getriebe<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>6.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>1</div><div>Übersetzung Getriebe Antrieb Zapfwelle<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>i_Zapfwelle<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>2.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>1</div><div>Wirkungsgrad Getriebe AntriebAbtrieb<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>eta_Getriebe<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>0.98<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>1</div><div>Wirkungsgrad Getriebe AntriebZapfwelle<span class=\"Apple-tab-span\" style=\"white-space:pre\">	</span>eta_Zapfwelle<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>0.95<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>1</div><div>Radius des Rades<span class=\"Apple-tab-span\" style=\"white-space:pre\">					</span>r_Rad<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>0.9<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>m</div><div>Masse des Rades<span class=\"Apple-tab-span\" style=\"white-space:pre\">					</span>m_Rad<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>80.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>kg</div><div>Masse des Chassis<span class=\"Apple-tab-span\" style=\"white-space:pre\">					</span>m_Chassis<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>1500.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>kg</div><div>Luftwiderstandsbeiwert Chassis<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>cW_Chassis<span class=\"Apple-tab-span\" style=\"white-space:pre\">		</span>0.95<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>1</div><div>Fläche Chassis Frontsicht<span class=\"Apple-tab-span\" style=\"white-space:pre\">				</span>A_Chassis<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>5.00<span class=\"Apple-tab-span\" style=\"white-space:pre\">			</span>m2</div><div><br></div><div><b>Beschreibung Simulationsmodell:</b></div><div>Das Motormodell stellt die antreibende Einheit des Traktors dar. Basierend auf einem Vergleich der Soll-Geschwindigkeit (velocity_Setpoint) und der Ist-Geschwindigkeit (vehicle_Velocity) wird von einem PID-Regler im Motor ein Motormoment erzeugt. Dieses Motormoment wiederum wird mit den Verlusten im Motor (abhängig von der Drehzahl) und dem resultierenden Moment aufgrund der Trägheit der Motorwelle verrechnet und gibt letztendlich ein Abtriebsmoment an den rotatorischen Connector des Motors.&nbsp;</div><div><br></div><div>Das Moment des Motors wird an das Getriebe weitergegeben, wo je nach Einstellung von Gang und Zapfwelle das Moment über die Parameter des Getriebes in ein Abtriebsmoment und Zapfwellenmoment umgewandelt wird.&nbsp;</div><div><br></div><div>Das Abtriebsmoment wiederum treibt das Rad des Traktors an, welches die rotatorische Größe (Moment und Winkel) in eine translatorische Größe (Kraft und Weg) übersetzt. Diese Kraft wird wiederum an das Chassis übergeben.</div><div><br></div><div>Das Chassis wiederum verrechnet die ankommende Kraft mit der Trägheit der eigenen Masse und der Kraft welche an der Anhängerkupplung anliegt. Im Falle, dass hier nichts angehängt ist, ist diese Kraft folglich = 0 N. Für den Fall, dass sich hier ein Arbeitsgerät befindet, wird der Weg den das Chassis zurücklegt weitergeleitet.</div><div><br></div><div>Im Arbeitsgerät wird die Wegänderung entsprechend in eine Beschleunigung umgerechnet, welche wieder für eine Kraft aufgrund der Trägheit des Arbeitsgerätes sorgt. Diese Kraft wird wieder über die Kupplung an die Zugmaschine gegeben und dort mit in die Beschleunigung des Traktors eingerechnet.</div><div><br></div><div><b>Wichtig:</b></div><div>Sofern an den Connectoren des Traktors keine entsprechenden Modelle angeschlossen sind, wird dort logischerweise auch keine Flussgröße übertragen. So kann die Kupplung des Traktors zwar Weg zurücklegen (Potenzialgröße), aber ohne entsprechendes Arbeitsgerät wird natürlich keine Kraft übertragen. Gleiches gilt für die Zapfwelle und die elektrischen Anschlüsse. Dort können zwar Rotationen stattfinden /Spannungen anliegen, aber ohne ein entsprechndes Modell, wird weder ein Moment, noch ein Strom übertragen.</div></body></html>", __OpenModelica_infoHeader = "<html><head></head><body>Das Modell \"Traktor\" ist die vereinfachte darstellung einer Landwirtschaftlichen Zugmaschine. Der Traktor kann in einem Modelica-Modell instanziiert werden. Der Traktor verfügt auf der Rückseite über Connectoren, welche das Anbinden von verschiedenen Arbeitsgeräten ermöglichen sollen.<div><br></div><div>Zusätzlich verfügt der Traktor über zwei Arten von Parametern:</div><div><ul><li>Während der Simulation veränderliche Parameter:</li><ul><li>Zielgeschwindigkeit (Real)</li><li>Gang (Boolscher Wert: Leerlauf = 0 &nbsp;oder Vorwärtsgang = 1)</li><li>Zapfwellen (Boolscher Wert: aktiv = 1 oder inaktiv = 0)</li><li>Batterieschalter (Boolscher Wert: On = 1 oder Off = 0)</li></ul></ul><br><ul><li>Während der Simulation unveränderliche Parameter</li><ul><li>Hierzu zählen die unterschiedlichen Parameter Teilmodelle, welche weiter unten näher erläutert werden.</li></ul></ul><br></div></body></html>"));
    end Traktor;
  end Zugmaschinen;

  package Modelle
    model Simulation
      Landwirtschaft.Zugmaschinen.Traktor traktor annotation(
        Placement(visible = true, transformation(origin = {44, 0}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 1, height = 11.11) annotation(
        Placement(visible = true, transformation(origin = {-65, 47}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep(startValue = false) annotation(
        Placement(visible = true, transformation(origin = {-64, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep1 annotation(
        Placement(visible = true, transformation(origin = {-64, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep2 annotation(
        Placement(visible = true, transformation(origin = {-64, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Landwirtschaft.Arbeitsgeraete.Masse masse annotation(
        Placement(visible = true, transformation(origin = {151, 1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    equation
      connect(booleanStep.y, traktor.Zapfwellensignal) annotation(
        Line(points = {{-52, -4}, {-26, -4}, {-26, 18}, {2, 18}}, color = {255, 0, 255}));
      connect(booleanStep1.y, traktor.Batterieschalter) annotation(
        Line(points = {{-52, -34}, {-16, -34}, {-16, 12}, {2, 12}}, color = {255, 0, 255}));
      connect(booleanStep2.y, traktor.Steuersingal_Input) annotation(
        Line(points = {{-52, -66}, {-10, -66}, {-10, 4}, {2, 4}}, color = {255, 0, 255}));
      connect(ramp.y, traktor.vel_Setpoint) annotation(
        Line(points = {{-46, 48}, {-20, 48}, {-20, 30}, {2, 30}}, color = {0, 0, 127}));
  connect(traktor.anhaengerkupplung, masse.translatoric) annotation(
        Line(points = {{82, 10}, {152, 10}, {152, 2}}));
    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation;

    model Testumgebung
    equation

      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "T", fontSize = 150)}));
    end Testumgebung;
    
    model Simulation_Heatloss
      Landwirtschaft.Zugmaschinen.Traktor traktor annotation(
        Placement(visible = true, transformation(origin = {-8, 12}, extent = {{-46, -46}, {46, 46}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 50, height = 10) annotation(
        Placement(visible = true, transformation(origin = {-117, 59}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep(startValue = false) annotation(
        Placement(visible = true, transformation(origin = {-116, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep1 annotation(
        Placement(visible = true, transformation(origin = {-116, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep2 annotation(
        Placement(visible = true, transformation(origin = {-116, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.LossyGear lossyGear(locked(fixed = true, start = true),lossTable = [0, 0.8, 0.8, 0, 0], mode(start = 1), startBackward(fixed = false), startForward(fixed = true, start = true), tau_etaPos(fixed = true, start = true), useHeatPort = true)  annotation(
        Placement(visible = true, transformation(origin = {96, 18}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Components.Convection convection annotation(
        Placement(visible = true, transformation(origin = {98, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature fixedTemperature(T = 25)  annotation(
        Placement(visible = true, transformation(origin = {138, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 20) annotation(
        Placement(visible = true, transformation(extent = {{28, -68}, {48, -48}}, rotation = 0)));
    Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe annotation(
        Placement(visible = true, transformation(origin = {160, 16}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
    equation
    connect(booleanStep.y, traktor.Zapfwellensignal) annotation(
        Line(points = {{-105, 8}, {-79, 8}, {-79, 30}, {-51, 30}}, color = {255, 0, 255}));
    connect(booleanStep1.y, traktor.Batterieschalter) annotation(
        Line(points = {{-105, -22}, {-69, -22}, {-69, 24}, {-51, 24}}, color = {255, 0, 255}));
    connect(booleanStep2.y, traktor.Steuersingal_Input) annotation(
        Line(points = {{-105, -54}, {-63, -54}, {-63, 16}, {-51, 16}}, color = {255, 0, 255}));
    connect(ramp.y, traktor.vel_Setpoint) annotation(
        Line(points = {{-98.3, 59}, {-72.3, 59}, {-72.3, 41}, {-50.3, 41}}, color = {0, 0, 127}));
    connect(traktor.zapfwelle, lossyGear.flange_a) annotation(
        Line(points = {{30.3333, 17.8267}, {74.3333, 17.8267}}));
    connect(const.y, convection.Gc) annotation(
        Line(points = {{49, -58}, {98, -58}, {98, -28}}, color = {0, 0, 127}));
    connect(fixedTemperature.port, convection.fluid) annotation(
        Line(points = {{148, -38}, {108, -38}}, color = {191, 0, 0}));
    connect(convection.solid, lossyGear.heatPort) annotation(
        Line(points = {{88, -38}, {88, -20}, {74, -20}, {74, -4}}, color = {191, 0, 0}));
    connect(lossyGear.flange_b, drehscheibe.D_Abtriebsseite) annotation(
        Line(points = {{118, 18}, {160, 18}, {160, 16}}));
    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_Heatloss;

    model Simulation_Kegelrad
      Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {67, 3}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Motor_mech motor_mech annotation(
        Placement(visible = true, transformation(origin = {-34, -16}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(motor_mech.abtrieb, kegelradgetriebe.K_Antriebsseite) annotation(
        Line(points = {{-9.76, -16.24}, {26.24, -16.24}, {26.24, -18.24}}));
      connect(step.y, motor_mech.velocity_setPoint) annotation(
        Line(points = {{-66, 30}, {-48, 30}, {-48, -10}}, color = {0, 0, 127}));
      connect(ramp.y, motor_mech.vehicle_velocity) annotation(
        Line(points = {{-67, -54}, {-48, -54}, {-48, -22}}, color = {0, 0, 127}));
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_Kegelrad;

    model Simulation_Kegelrad_Drehscheibe
      Landwirtschaft.Hilfsmodelle.Motor_mech motor_mech annotation(
        Placement(visible = true, transformation(origin = {-34, -16}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {58, -2}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step1(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step2(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe1 annotation(
        Placement(visible = true, transformation(origin = {36, 62}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
      // Gleichungen
    equation
//zahnrad.Z_Antriebsseite.alpha = kegelradgetriebe.K_Abtriebsseite.alpha;
      connect(step.y, motor_mech.velocity_setPoint) annotation(
        Line(points = {{-66, 30}, {-48, 30}, {-48, -10}}, color = {0, 0, 127}));
      connect(motor_mech.abtrieb, kegelradgetriebe.K_Antriebsseite) annotation(
        Line(points = {{-10, -16}, {6, -16}, {6, -17}, {28, -17}}));
      connect(step2.y, motor_mech.vehicle_velocity) annotation(
        Line(points = {{-66, -58}, {-48, -58}, {-48, -22}}, color = {0, 0, 127}));
      connect(kegelradgetriebe.K_Abtriebsseite, drehscheibe1.D_Abtriebsseite) annotation(
        Line(points = {{66, 20}, {36, 20}, {36, 62}}));
    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_Kegelrad_Drehscheibe;

    model Simulation_Kegelrad_Drehscheibe_Zahnrad
      Landwirtschaft.Hilfsmodelle.Motor_mech motor_mech annotation(
        Placement(visible = true, transformation(origin = {-34, -16}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {58, -2}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step1(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step2(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Zahnrad zahnrad annotation(
        Placement(visible = true, transformation(origin = {65, 55}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe annotation(
        Placement(visible = true, transformation(origin = {-13, 65}, extent = {{-33, -33}, {33, 33}}, rotation = 0)));
      // Gleichungen
    equation
//zahnrad.Z_Antriebsseite.alpha = kegelradgetriebe.K_Abtriebsseite.alpha;
      connect(step.y, motor_mech.velocity_setPoint) annotation(
        Line(points = {{-66, 30}, {-48, 30}, {-48, -10}}, color = {0, 0, 127}));
      connect(motor_mech.abtrieb, kegelradgetriebe.K_Antriebsseite) annotation(
        Line(points = {{-10, -16}, {6, -16}, {6, -17}, {28, -17}}));
      connect(step2.y, motor_mech.vehicle_velocity) annotation(
        Line(points = {{-66, -58}, {-48, -58}, {-48, -22}}, color = {0, 0, 127}));
      connect(kegelradgetriebe.K_Abtriebsseite, zahnrad.Z_Antriebsseite) annotation(
        Line(points = {{66, 20}, {36, 20}, {36, 56}, {58, 56}}));
      connect(zahnrad.Z_Abtriebsseite, drehscheibe.D_Abtriebsseite) annotation(
        Line(points = {{72, 56}, {90, 56}, {90, 94}, {-12, 94}, {-12, 66}}));
    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_Kegelrad_Drehscheibe_Zahnrad;

    model Simulation_Kegelrad_Drehscheibe12_Zahnrad
      Landwirtschaft.Hilfsmodelle.Motor_mech motor_mech annotation(
        Placement(visible = true, transformation(origin = {-34, -16}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {58, -2}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
      Modelica.Blocks.Sources.Step step2(height = 10) annotation(
        Placement(visible = true, transformation(origin = {-78, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Zahnrad zahnrad annotation(
        Placement(visible = true, transformation(origin = {65, 55}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe2 annotation(
        Placement(visible = true, transformation(origin = {65, -65}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
      Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe1 annotation(
        Placement(visible = true, transformation(origin = {-21, 69}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
      // Gleichungen
    equation
//zahnrad.Z_Antriebsseite.alpha = kegelradgetriebe.K_Abtriebsseite.alpha;
      connect(step.y, motor_mech.velocity_setPoint) annotation(
        Line(points = {{-66, 30}, {-48, 30}, {-48, -10}}, color = {0, 0, 127}));
      connect(motor_mech.abtrieb, kegelradgetriebe.K_Antriebsseite) annotation(
        Line(points = {{-10, -16}, {6, -16}, {6, -17}, {28, -17}}));
      connect(step2.y, motor_mech.vehicle_velocity) annotation(
        Line(points = {{-66, -58}, {-48, -58}, {-48, -22}}, color = {0, 0, 127}));
      connect(kegelradgetriebe.K_Abtriebsseite, zahnrad.Z_Antriebsseite) annotation(
        Line(points = {{66, 20}, {36, 20}, {36, 56}, {58, 56}}));
      connect(kegelradgetriebe.K_Abtriebsseite, drehscheibe2.D_Abtriebsseite) annotation(
        Line(points = {{66, 20}, {90, 20}, {90, -64}, {66, -64}}));
      connect(zahnrad.Z_Abtriebsseite, drehscheibe1.D_Abtriebsseite) annotation(
        Line(points = {{72, 56}, {90, 56}, {90, 90}, {-20, 90}, {-20, 70}}));
    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_Kegelrad_Drehscheibe12_Zahnrad;

    model Simulation_Traktor_Kegelrad
      Landwirtschaft.Zugmaschinen.Traktor traktor annotation(
        Placement(visible = true, transformation(origin = {-3, 3}, extent = {{-37, -37}, {37, 37}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep annotation(
        Placement(visible = true, transformation(origin = {-76, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep1 annotation(
        Placement(visible = true, transformation(origin = {-76, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.BooleanStep booleanStep2 annotation(
        Placement(visible = true, transformation(origin = {-76, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 15) annotation(
        Placement(visible = true, transformation(origin = {-80, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {79, 19}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Zahnrad zahnrad annotation(
        Placement(visible = true, transformation(origin = {-22, 74}, extent = {{-28, -28}, {28, 28}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe1(J_D = 200)  annotation(
        Placement(visible = true, transformation(origin = {18, 86}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe2(J_D = 200)  annotation(
        Placement(visible = true, transformation(origin = {46, 84}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  Modelica.Blocks.Sources.Step step(height = 11.11, offset = 0)  annotation(
        Placement(visible = true, transformation(origin = {-80, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(booleanStep.y, traktor.Zapfwellensignal) annotation(
        Line(points = {{-64, -14}, {-58, -14}, {-58, 18}, {-38, 18}}, color = {255, 0, 255}));
      connect(booleanStep1.y, traktor.Batterieschalter) annotation(
        Line(points = {{-64, -44}, {-54, -44}, {-54, 12}, {-38, 12}}, color = {255, 0, 255}));
      connect(booleanStep2.y, traktor.Steuersingal_Input) annotation(
        Line(points = {{-64, -74}, {-48, -74}, {-48, 6}, {-38, 6}}, color = {255, 0, 255}));
  connect(traktor.zapfwelle, kegelradgetriebe.K_Antriebsseite) annotation(
        Line(points = {{28, 8}, {40, 8}, {40, 7}, {54, 7}}));
  connect(kegelradgetriebe.K_Abtriebsseite, zahnrad.Z_Antriebsseite) annotation(
        Line(points = {{85, 38}, {-28, 38}, {-28, 73}}));
  connect(kegelradgetriebe.K_Abtriebsseite, drehscheibe2.D_Abtriebsseite) annotation(
        Line(points = {{85, 38}, {85, 84}, {46, 84}}));
  connect(zahnrad.Z_Abtriebsseite, drehscheibe1.D_Abtriebsseite) annotation(
        Line(points = {{-16, 73}, {18, 73}, {18, 86}}));
  connect(step.y, traktor.vel_Setpoint) annotation(
        Line(points = {{-68, 56}, {-56, 56}, {-56, 28}, {-38, 28}}, color = {0, 0, 127}));
    annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));end Simulation_Traktor_Kegelrad;
    
    model Simulation_Tank
    Landwirtschaft.Hilfsmodelle.Tank tank(V(start = 10))  annotation(
        Placement(visible = true, transformation(origin = {1, -1}, extent = {{-59, -59}, {59, 59}}, rotation = 0)));
    equation

    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_Tank;
    
    model Simulation_FirstRunningPrototyp
    Landwirtschaft.Zugmaschinen.Traktor traktor annotation(
        Placement(visible = true, transformation(origin = {-14, 2}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
  Modelica.Blocks.Sources.Step step(height = 11)  annotation(
        Placement(visible = true, transformation(origin = {-82, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep annotation(
        Placement(visible = true, transformation(origin = {-84, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep1 annotation(
        Placement(visible = true, transformation(origin = {-84, -24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanStep booleanStep2 annotation(
        Placement(visible = true, transformation(origin = {-82, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Landwirtschaft.Arbeitsgeraete.Anhaenger anhaenger annotation(
        Placement(visible = true, transformation(origin = {105, -3}, extent = {{-51, -51}, {51, 51}}, rotation = 0)));
  Landwirtschaft.Arbeitsgeraete.Anhaenger_Rad anhaenger_Rad annotation(
        Placement(visible = true, transformation(origin = {122, -56}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Kegelradgetriebe kegelradgetriebe annotation(
        Placement(visible = true, transformation(origin = {75, 57}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Zahnrad zahnrad annotation(
        Placement(visible = true, transformation(origin = {130, 52}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe1 annotation(
        Placement(visible = true, transformation(origin = {131, 85}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Drehscheibe drehscheibe2 annotation(
        Placement(visible = true, transformation(origin = {171, 59}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
  Landwirtschaft.Hilfsmodelle.Tank tank annotation(
        Placement(visible = true, transformation(origin = {1, 65}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    equation
      connect(step.y, traktor.vel_Setpoint) annotation(
        Line(points = {{-70, 48}, {-62, 48}, {-62, 28}, {-52, 28}}, color = {0, 0, 127}));
  connect(booleanStep.y, traktor.Zapfwellensignal) annotation(
        Line(points = {{-72, 10}, {-64, 10}, {-64, 18}, {-52, 18}}, color = {255, 0, 255}));
  connect(booleanStep1.y, traktor.Batterieschalter) annotation(
        Line(points = {{-72, -24}, {-66, -24}, {-66, 12}, {-52, 12}}, color = {255, 0, 255}));
  connect(booleanStep2.y, traktor.Steuersingal_Input) annotation(
        Line(points = {{-70, -58}, {-60, -58}, {-60, 6}, {-52, 6}}, color = {255, 0, 255}));
  connect(traktor.anhaengerkupplung, anhaenger.anhaenger_kupplung) annotation(
        Line(points = {{20, 12}, {34, 12}, {34, 2}, {58, 2}}));
  connect(anhaenger_Rad.radlager, anhaenger.translatoric_radlager) annotation(
        Line(points = {{90, -52}, {90, -26}, {128, -26}, {128, -12}}));
  connect(traktor.zapfwelle, kegelradgetriebe.K_Antriebsseite) annotation(
        Line(points = {{20, 8}, {36, 8}, {36, 48}, {56, 48}}));
  connect(kegelradgetriebe.K_Abtriebsseite, zahnrad.Z_Antriebsseite) annotation(
        Line(points = {{80, 70}, {104.5, 70}, {104.5, 52}, {125, 52}}));
  connect(kegelradgetriebe.K_Abtriebsseite, drehscheibe1.D_Abtriebsseite) annotation(
        Line(points = {{80, 70}, {78, 70}, {78, 86}, {132, 86}}));
  connect(zahnrad.Z_Abtriebsseite, drehscheibe2.D_Abtriebsseite) annotation(
        Line(points = {{134, 52}, {172, 52}, {172, 60}}));
  connect(tank.translatoric_tank, traktor.anhaengerkupplung) annotation(
        Line(points = {{1, 65}, {28, 65}, {28, 12}, {20, 12}}));
    protected
      annotation(
        Icon(graphics = {Text(extent = {{-100, 100}, {100, -100}}, textString = "S", fontSize = 150)}));
    end Simulation_FirstRunningPrototyp;
  end Modelle;
  annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Bitmap(origin = {-4, 8}, extent = {{88, -158}, {-88, 158}}, fileName = "modelica://Landwirtschaft/Bib_Image.png")}));
end Landwirtschaft;
