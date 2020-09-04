%coordinate setup for CO relaxation 
%define C(i), CA(i), O(i), N(i+1), HA(i), NH(i+1), CA(i+1),HA(i+1)

function [cphiCcsax, cthetaCcsax, cphiCcsay, cthetaCcsay, cphiNcsax, cthetaNcsax, cphiNcsay, cthetaNcsay, cphiNH, cthetaNH, cphiCCa, cthetaCCa, cphiCH, cthetaCH, cphiCN, cthetaCN, cphiNCA, cthetaNCA] = coords_global_peptide_planes(H, N, CA, C, O)
	CO=O-C;
	CN=N-C;
	 CCA=CA-C;
	% CHA=HA-C;
	CH=H-C;
	HN=H-N;
	% HACA=HA-CA;
	NCA=N-CA;
	% HAp1CAp1=HAp1-CAp1;
	for a=1:length(CO)
		CO(a,:)=CO(a,:)/norm(CO(a,:));
		CN(a,:)=CN(a,:)/norm(CN(a,:));
		CCA(a,:)=CCA(a,:)/norm(CCA(a,:));
	%     CHA(a,:)=CHA(a,:)/norm(CHA(a,:));
		CH(a,:)=CH(a,:)/norm(CH(a,:));
		HN(a,:)=HN(a,:)/norm(HN(a,:));
	%     HACA(a,:)=HACA(a,:)/norm(HACA(a,:));
		NCA(a,:)=NCA(a,:)/norm(NCA(a,:));
	%     HAp1CAp1(a,:)=HAp1CAp1(a,:)/norm(HAp1CAp1(a,:));
	end
		
	%coordinates for CSA relaxation
	for a=1:length(HN);
		%for i=1:1
		Sxx(a,1:3)=rotateIt(cross(CO(a,:),CN(a,:)),36,CN(a,:));
		cphiCcsax(a)=atan((Sxx(a,2))./(Sxx(a,1)));
		cthetaCcsax(a)=acos(Sxx(a,3)./(sqrt(Sxx(a,1)^2+Sxx(a,2)^2+Sxx(a,3)^2)));
		Syy(a,1:3)=rotateIt(cross(CO(a,:),CN(a,:)),-54,CN(a,:));
		cphiCcsay(a)=atan((Syy(a,2))./(Syy(a,1)));
		cthetaCcsay(a)=acos(Syy(a,3)/sqrt(Syy(a,1)^2+Syy(a,2)^2+Syy(a,3)^2));



		%H(i+1)N(i+1) coordinates as well as sXX for 15N CSA
		% S11(a,1:3)=rodrigues_rot(HN(a,:),cross(HN(a,:),CN(a,:)),deg2rad(-17));
		% S33(a,1:3)=rodrigues_rot(HN(a,:),cross(HN(a,:),CN(a,:)),sym(pi/2)-deg2rad(-17));
		S11(a,1:3)=rotateIt(cross(CO(a,:),CN(a,:)),-(58+17),CN(a,:));
		S33(a,1:3)=rotateIt(cross(CO(a,:),CN(a,:)),-(58+17-90),CN(a,:));
		% S11(a,1:3)=rotateIt(cross(CO(a,:),CN(a,:)),-17,HN(a,:));
		%Sn(a,1:3)=HN(a,:);
		% Sn(a,1:3)=rodrigues_rot(HN(a,:),cross(CO(a,:),CN(a,:)),deg2rad(-17));
		% Sn2(a,1:3)=rotateIt(cross(CO(a,:),CN(a,:)),-90,HN(a,:));
		cphiNcsax(a)=atan(S11(a,2)/(S11(a,1)));
		cthetaNcsax(a)=acos(S11(a,3)/sqrt(S11(a,1)^2+S11(a,2)^2+S11(a,3)^2));
		% Sn3(a,:)=cross(Sn(a,:),Sn2(a,:));
		% S22(a,:)=rotateIt(S33(a,:),90,S11(a,:));
		%the choice of the two axially symmetric components is arbitrary, chose S33
		%instead of S22 to avoid rounding errors in Matlab associated with
		%approximations of rotations; update the relaxation rate accordingly
		%(exchange x and y)
		cphiNcsay(a)=atan(S33(a,2)/S33(a,1));
		cthetaNcsay(a)=acos(S33(a,3)/sqrt(S33(a,1)^2+S33(a,2)^2+S33(a,3)^2));
		% S2CCSAxy=GAF_ord_param(5,5,5,cphiCcsax(a),cthetaCcsax(a),cphiCcsay(a),cthetaCcsay(a))
		% if (S2CCSAxy >= -0.5) && (S2CCSAxy <= 0)
		%     Sn3(a,:)=cross(Sn(a,:),Sn2(a,:));
		% else
		%     Sn3(a,:)=cross(Sn2(a,:),Sn(a,:));
		% end

		cphiNH(a)=atan(HN(a,2)/HN(a,1));
		cthetaNH(a)=acos(HN(a,3)/sqrt(HN(a,1)^2+HN(a,2)^2+HN(a,3)^2));

		%C'(i)CA(i)
		cphiCCa(a)=atan(CCA(a,2)/CCA(a,1));
		cthetaCCa(a)=acos(CCA(a,3)/sqrt(CCA(a,1)^2+CCA(a,2)^2+CCA(a,3)^2)); % check me

		%CA(i)HA(i)
		% cphiCAHA(a)=atan(HACA(a,2)/HACA(a,1));
		% cthetaCAHA(a)=acos(HACA(a,3)/sqrt(HACA(a,1)^2+HACA(a,2)^2+HACA(a,3)^2));

		%C'(i)H(i+1)
		cphiCH(a)=atan(CH(a,2)/CH(a,1));
		cthetaCH(a)=acos(CH(a,3)/sqrt(CH(a,1)^2+CH(a,2)^2+CH(a,3)^2));

		%C(i)HA(i)
		% cphiCHA(a)=atan(CHA(a,2)/CHA(a,1));
		% cthetaCHA(a)=acos(CHA(a,3)/sqrt(CHA(a,1)^2+CHA(a,2)^2+CHA(a,3)^2));

		%C'(i)N(i+1)
		cphiCN(a)=atan(CN(a,2)/CN(a,1));
		cthetaCN(a)=acos(CN(a,3)/sqrt(CN(a,1)^2+CN(a,2)^2+CN(a,3)^2));

		%CA(i+1)N(i+1)
		cphiNCA(a)=atan(NCA(a,2)/NCA(a,1));
		cthetaNCA(a)=acos(NCA(a,3)/sqrt(NCA(a,1)^2+NCA(a,2)^2+NCA(a,3)^2));


		
		%%cphiCaN(a)=atan(NCA(a,2)/NCA(a,1));
		%%cthetaCan(a)=acos(NCA(a,3)/sqrt(NCA(a,1)^2+NCA(a,2)^2+NCA(a,3)^2));

	end
	
	
	
end

	%save COCSA.mat res cphi ctheta CSiso R1C800
	% Sxx(i,1:3)=rotateIt(cross(CO(i,:),CN(i,:)),36,CN(i,:));
	% cphix2(i)=acos((Sxx(i,1))./(sqrt((Sxx(i,1)).^2+(Sxx(i,2)).^2)));
	% cthetax2(i)=acos((Sxx(i,3))./(sqrt((Sxx(i,1)).^2+(Syy(i,2)).^2+(Syy(i,3)).^2)));
	% Syy(i,1:3)=rotateIt(cross(CO(i,:),CN(i,:)),55.4,CCA(i,:));
	% cphiy2(i)=acos((Syy(i,1))./(sqrt((Syy(i,1)).^2+(Syy(i,2)).^2)));
	% cthetay2(i)=acos((Syy(i,3))./(sqrt((Syy(i,1)).^2+(Syy(i,2)).^2+(Syy(i,3)).^2)));
	% cphix3(i)=atan2(Sxx(i,2),Sxx(i,1))
