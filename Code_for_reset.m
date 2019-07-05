clear all
close all
clc

scale_ed=ones(1,6); scale_ped=ones(1,2);
scale_md=ones(1,6); scale_pmd=ones(1,2);
scale_ch=ones(1,6); scale_om=ones(1,2);
scale_te=ones(1,6); scale_mo=ones(1,2);
scale_ch_max=1;     scale_te_max=1;


theta_ed=zeros(1,6); theta_ped=zeros(1,2);
theta_md=zeros(1,6); theta_pmd=zeros(1,2);
theta_ch=zeros(1,6); theta_om=zeros(1,2);
theta_te=zeros(1,6); theta_mo=zeros(1,2);
theta_opt_ch=0;      theta_opt_te=0;

phi_ed=zeros(1,6); phi_ped=zeros(1,2);
phi_md=zeros(1,6); phi_pmd=zeros(1,2);
phi_ch=zeros(1,6); phi_om=zeros(1,2);
phi_te=zeros(1,6); phi_mo=zeros(1,2);
phi_opt_ch=0;      phi_opt_te=0;

handedness_ch=ones(1,6);
orientation_te=ones(1,6);
handedness_ch_max=1;
orientation_te_max=1;

fileID = fopen('resetting the modules.py','w');
ProjectName='Materiatronics_modules';
fprintf(fileID, 'import ScriptEnv\n');
fprintf(fileID, 'ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")\n');
fprintf(fileID, 'oDesktop.RestoreWindow()\n');
fprintf(fileID, 'oProject = oDesktop.SetActiveProject("%s")\n',ProjectName);
fprintf(fileID, 'oDesign = oProject.SetActiveDesign("HFSSDesign1")\n');
fprintf(fileID, 'oEditor = oDesign.SetActiveEditor("3D Modeler")\n');

% ee reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Electric_dipole_%s:Scale:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ed(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ed(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ed(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Electric_dipole_%s:Rotate:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_ed(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Electric_dipole_%s:Rotate:2"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_ed(ic)));
end;

% ee non-reciprocal
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_real:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ped(1)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ped(1)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ped(1)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_real:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_ped(1))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_real:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_ped(1))));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_imaginary:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ped(2)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ped(2)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ped(2)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_imaginary:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_ped(2))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_electric_dipole_imaginary:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_ped(2))));


% mm reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Scale:1"],["NAME:ChangedProps",',num2str(ic*2+2));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_md(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_md(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_md(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:1"],',num2str(ic*2+2));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_md(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:2"],',num2str(ic*2+2));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_md(ic)));
end;
 
% mm non-reciprocal
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_real:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_pmd(1)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_pmd(1)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_pmd(1)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_real:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_pmd(1))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_real:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_pmd(1))));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_imaginary:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_pmd(2)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_pmd(2)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_pmd(2)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_imaginary:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_pmd(2))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Precessing_magnetic_dipole_imaginary:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_pmd(2))));


% em reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:Scale:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ch(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ch(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ch(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:Rotate:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_ch(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:Rotate:2"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_ch(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve%s:CreateEquationCurve:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Y(_t)","Value:=","%slambda/25*sin(_t)"]]]])\n',[num2str(handedness_ch(ic)) '*']);
end;

for ic=1:2
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Scale:1"],["NAME:ChangedProps",',num2str(ic*2+22));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_om(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_om(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_om(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:3"],',num2str(ic*2+22));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_om(ic))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Cylinder%s:Rotate:4"],',num2str(ic*2+22));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_om(ic))));
end;


% em non-reciprocal
for ic=1:6
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:Scale:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_te(ic)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_te(ic)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_te(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:Rotate:1"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_te(ic)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:Rotate:2"],',num2str(ic));
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_te(ic)));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptel%s:CreateCylinder:1"],["NAME:ChangedProps",',num2str(ic));
fprintf(fileID, '["NAME:Center Position","X:=","0","Y:=","%slambda/25","Z:=", "-lambda/5/2"]]]])\n',[num2str(orientation_te(ic)) '*']);
end;
 
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_real:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_mo(1)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_mo(1)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_mo(1)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_real:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_mo(1))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_real:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_mo(1))));

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_imag:Scale:2"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_mo(2)));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_mo(2)));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_mo(2)));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_imag:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(theta_mo(2))));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","Moving_imag:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(round(phi_mo(2))));

% orientation for maximum chirality and Tellegen effects
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:Scale:1"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_ch_max));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_ch_max));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_ch_max));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_opt_ch));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_opt_ch));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","EquationCurve7:CreateEquationCurve:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Y(_t)","Value:=","%slambda/25*sin(_t)"]]]])\n',[num2str(handedness_ch_max) '*']);

fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:Scale:1"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Scale of X","Value:=", "%s"],',num2str(scale_te_max));
fprintf(fileID, '["NAME:Scale of Y","Value:=", "%s"],',num2str(scale_te_max));
fprintf(fileID, '["NAME:Scale of Z","Value:=", "%s"]]]])\n',num2str(scale_te_max));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:Rotate:1"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(theta_opt_te));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:Rotate:2"],');
fprintf(fileID, '["NAME:ChangedProps",["NAME:Angle",	"Value:="," %sdeg"]]]])\n',num2str(phi_opt_te));
fprintf(fileID, 'oEditor.ChangeProperty(["NAME:AllTabs",["NAME:Geometry3DCmdTab",["NAME:PropServers","striptelmax:CreateCylinder:1"],["NAME:ChangedProps",');
fprintf(fileID, '["NAME:Center Position","X:=","0","Y:=","%slambda/25","Z:=", "-lambda/5/2"]]]])\n',[num2str(orientation_te_max) '*']);

fclose(fileID)

