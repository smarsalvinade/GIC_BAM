% Property: Observatori de l'Ebre. Authors: Victoria Canillas Pérez and Santiago Marsal Vinadé.
% License: GNU General Public License 3 (GPL-3.0) (http://opensource.org/licenses/GPL-3.0)

% Subprogram to construct a DC-equivalent network from sheets 'Substations' and 'Lines' of the defined Excel file. It is called by the GIC_BAM main program.

eps = 1e-6; % Small value assigned to zero earthing and winding resistances [ohm].
% Determine the vector of sorted bus voltage levels in 'Substations' sheet:
VLTF = [HVBusTF, LVBusTF]; % High and low voltage levels of each transformer (nTF x 2) [kVAC]
if any(isnan(HVBusTF)); errordlg('All transformers should have a HV value (column 11 of ''Substations'' sheet.', 'File error'); end % Control
VLSubs = NaN(max(nTFSubs*2), m); % Will contain the sorted voltage levels at each substation (L x m) [kVAC].
for Subs = 1: m % Loop through substations to set the sorted voltage levels at each substation
    RowsCurSubs = RowsSubs(Subs) + (0: nTFSubs(Subs) - 1); % Rows in 'Substations' sheet corresponding to the current substation (excluding header).
    VLTFSubs = unique(VLTF(RowsCurSubs, :)); VLTFSubs = VLTFSubs(~isnan(VLTFSubs)); VLTFSubs = VLTFSubs(:); % Voltage levels found at current substation [kVAC].
    VLSubs(1: length(VLTFSubs), Subs) = flipud(VLTFSubs);
end
l = sum(~isnan(VLSubs)); % Nº of voltage levels providing a path for GIC in each substation (1 x m).
L = max(l); % Maximum nº of voltage levels in a single substation. This is the allocated number of buses in each substation of the model.
VLSubs = VLSubs(1: L, :); % Sorted (descending order) voltage levels at each substation (L x m) [kVAC].
nNodCon = L*(L + 1)/2; % Nº of possible connections between nodes of the same substation.
ntot = n + nNodCon*m; % Total nº of branches (real TLs + subs. internal branches) in the DC-equivalent network.
nNod = m*(L + 1); % Total nº of nodes in the equivalent network.
nBBCon = L*(L - 1)/2; % nº of possible connections between buses of the same substation.

% Compute LinesEq matrix containing information on the lines of the DC-equivalent network:
LinesEq = NaN(ntot, 8); LinesEq(1: n, [3, 4, 6, 7, 8]) = [RphLin, nrepLin, VTLin, OriSubsLin, DesSubsLin]; % (ntot x 8)
for tlin = 1: n % Loop through (real) transmission lines to set their origin and destination nodes in the equivalent network.
    VLin = VTLin(tlin); % Line voltage [kVAC].
    SubsLineOri = OriSubsLin(tlin); SubsLineDes = DesSubsLin(tlin); % Nº of the substation where the current line starts and nº of the substation where the current line ends.
    if ~any(VLSubs(:, SubsLineOri) == VLin); errordlg(['Voltage of line nº ', num2str(tlin), ' (according to ''Lines'' sheet) not found among voltages of its origin substation (according to ''Substations'' sheet).'], 'File error'); end % Control
    if ~any(VLSubs(:, SubsLineDes) == VLin); errordlg(['Voltage of line nº ', num2str(tlin), ' (according to ''Lines'' sheet) not found among voltages of its destination substation (according to ''Substations'' sheet).'], 'File error'); end % Control
    LinesEq(tlin, 1: 2) = [m*find(VLin == VLSubs(:, SubsLineOri)) + SubsLineOri, m*find(VLin == VLSubs(:, SubsLineDes)) + SubsLineDes]; % Origin and destination nodes of transmission lines.
end

% Generate new nodes in NodesEq array, namely, L + 1 nodes per real substation:
% 1st node corresponds to substation neutral point (NP).
% 2nd node corresponds to highest AC voltage level bus (b1) found in substation.
% (l(Subs) + 1)-th node corresponds to lowest AC voltage level bus (bl) found in substation.
% (l(Subs) + 2)-th node corresponds to the first "vitrual" bus.
% (L + 1)-th node corresponds to the last "virtual" bus.
NodesEq = NaN(nNod, 7); % Will contain numerical info about the nodes of the equivalent network. (nNod x 7)
NodesEqTxt = cell(nNod, 7); % Will contain text info about the nodes of the equivalent network. (nNod x 7)
RgSubs(RgSubs == 0) = eps;
transVLSubs = VLSubs';
NodesEq(:, [2: 5, 7]) = [(1: nNod)', repmat([LatLonSubs, RgSubs], L + 1, 1), [NaN(m, 1); transVLSubs(:)]]; NodesEq(m + 1: end, 5) = NaN; % Node nº, Latitude [º], Longitude [º], earthing resistence [ohm] and voltage levels [kVAC].
NodesEqTxt(1: m, 1) = cellstr(repmat("np_Sub", m, 1) + NodesEq(1: m, 2)); % Names of the nodes corresponding to the neutral point of each substation (np_Subk).
for i = 1: L
    NodesEqTxt(i*m + 1: (i + 1)*m, 1) = cellstr(repmat("b", m, 1) + i + repmat("_Sub", m, 1) + NodesEq(1: m, 2)); % Names of the nodes corresponding to the buses (bi_Subk).
end
RphW1STF(RphW1STF == 0) = eps; RphW1Sinv = 1./RphW1STF; RphW1Sinv(isnan(RphW1Sinv)) = 0; % Resistance of HV (GSU, GY-GY) or series (A) winding of each transformer. (nTF x 1) [ohm/ph]
RphW2CTF(RphW2CTF == 0) = eps; RphW2Cinv = 1./RphW2CTF; RphW2Cinv(isnan(RphW2Cinv)) = 0; % Resistance of LV (GSU, GY-GY) or common (A) winding of each transformer. (nTF x 1) [ohm/ph]
IsATall = strcmp('A', TFtypeTF); % 1 for autotransformers; 0 otherwise. (nTF x 1)
IsTeeall = strcmp('T', TFtypeTF); % 1 for Tee substations; 0 otherwise. (nTF x 1)
for Subs = 1: m
    RowsCurSubs = RowsSubs(Subs) + (0: nTFSubs(Subs) - 1); % Rows in 'Substations' sheet corresponding to the current substation (excluding header).
    Reqbn = NaN(L, 1); % Will contain the equivalent resistance of bus-neutral connections (primary, secondary and common winding resistances of transformers). [ohm] (L x 1)
    ReqSer = NaN(nBBCon, 1); % Will contain the equivalent resistance of bus-bus connections within substation (series winding resistances of autotransformers). [ohm] (nBBCon x 1)
    HVBus = HVBusTF(RowsCurSubs); % HV bus of the transformers in current substation [kVAC]. Cannot be blank. (nTFSubs(Subs) x 1)
    LVBus = LVBusTF(RowsCurSubs); % LV bus of the transformers in current substation [kVAC]. Can be blank (e.g. for delta windings, which are not earthed, and tee's). (nTFSubs(Subs) x 1)
    IsATSubs = IsATall(RowsCurSubs); if any(isnan(HVBus(IsATSubs) + LVBus(IsATSubs))); errordlg(['Bus voltages cannot be empty for an autotransformer. Revise columns 11 and 12 of substation ', num2str(Subs), ' in ''Substations'' sheet.'], 'File error'); end  % 1 for autotransformers; 0 otherwise. (nTFSubs(Subs) x 1)
    IsTeeSubs = IsTeeall(RowsCurSubs);  % 1 for Tee's; 0 otherwise. (nTFSubs(Subs) x 1)
    RphW1SinvSubs = RphW1Sinv(RowsCurSubs); RphW2CinvSubs = RphW2Cinv(RowsCurSubs); % Reciprocal resistance of transformern windings. [ohm^-1] (nTFSubs(Subs) x 1)
    CountConNod = 0;
    for i = 0: L % Loop through the L + 1 nodes of each substation.
        CountConNod = CountConNod(end) + (1: L - i); % Nº of connections between node i and the other nodes of the same substation.
        LinesEq(n + nNodCon*(Subs - 1) + CountConNod, [1, 2, 7, 8]) = [repmat(i*m + Subs, L - i, 1), m*(i + 1: L)' + Subs, repmat(Subs, L - i, 2)]; % Origin node, destination node, origin substation and destination substation of the new branches joining the L + 1 nodes of a given substation.
        TxtRefTrafos = "";
        for tf = 1: nTFSubs(Subs)
            if i <= l(Subs) && (i <= 1 || any(VLSubs(i, Subs) == [HVBus(tf), LVBus(tf)]))
                TxtRefTrafos = TxtRefTrafos + TxtRefTF(RowsSubs(Subs) - 1 + tf) + ", "; % Trafo reference text corresponding to the substation associated with each node (to track original transformers).
            end
        end
        NodesEqTxt{i*m + Subs, 6} = TxtRefTrafos{1}(1: end - 2); % Write transformer reference text in NodesEq array.
        if i > 0 && i <= l(Subs)
            Reqbn(i) = 1/(sum(RphW1SinvSubs(HVBus == VLSubs(i, Subs) & ~IsATSubs)) + sum(RphW2CinvSubs(LVBus == VLSubs(i, Subs)))); % Equivalent resistance of bus-neutral connections (primary, secondary and common winding resistances of transformers). [ohm] (L x 1)
            if any(IsTeeSubs) % Tee substations:
                if nTFSubs(Subs) > 1; errordlg(['Tee substations are expected to consist of a single row in ''Substations'' sheet, i.e., no other transformers are allowed. Revise substation nº ', num2str(Subs), '.'], 'File error'); end % Control
                if any(~isnan([RgSubs(Subs), LVBus])); errordlg(['LV bus (column 12) and Rg (column 5) should be empty cells in a Tee substation. Revise substation nº ', num2str(Subs), ' in ''Substations'' sheet.'], 'File error'); end % Control
            end
            LinesEq(n + nNodCon*(Subs - 1) + i, 3) = Reqbn(i); % Equivalent bus-neutral resistance of the branches joining the L buses with the np at a given substation. [ohm] (L x 1)
            for j = i + 1: l(Subs)
                BBCon = L*(i - 1) - i*(i + 1)/2 + j; % Bus-bus connection.
                ReqSer(BBCon) = 1/(sum(RphW1SinvSubs(HVBus == VLSubs(i, Subs) & LVBus == VLSubs(j, Subs) & IsATSubs))); % Equivalent resistance of bus-bus connections within substation (series winding resistances of autotransformers). [ohm] (nBBCon x 1)
                LinesEq(n + nNodCon*(Subs - 1) + L + BBCon, 3) = ReqSer(BBCon); % Equivalent bus-bus series resistance of the branches joining the L buses of a given substation. [ohm] (nBBCon x 1)
            end
        end
    end
    if all(~isfinite(LinesEq(n + nNodCon*(Subs - 1) + (1: L), 3))) % If there's no way to ground in substation, allocate a finite grounding resistance while keeping infinite bus-neutral resistances.
        if ~isfinite(NodesEq(Subs, 5)); NodesEq(Subs, 5) = 0.10012; end
    end
end
LinesEq(n + 1: end, 4) = 1; LinesEq(:, 5) = LinesEq(:, 3)/3./[nrepLin; ones(ntot - n, 1)]; LinesEq(~isfinite(LinesEq)) = NaN;

% Write NodesEq and LinesEq arrays in Excel file if required:
if WriteEqNetXLS
    NodesEqTxt(:, [2: 5, 7]) = num2cell(NodesEq(:, [2: 5, 7])); NodesEqTxt = [{'Node', 'Nº', 'Latitude (°)', 'Longitude (°)', 'Rg (ohm)', 'Trafo ref.', 'Bus voltage (kVAC)'}; NodesEqTxt]; xlswrite(XLSInpFileName, NodesEqTxt, 'NodesEq');
    LinesEqTxt = [{'Origin node #', 'Destination node #', 'Resistance (ohm/ph)', 'Nº lines', 'Total line resistance (ohm)', 'Line voltage (kVAC)', 'Origin subs #', 'Destination subs #'}; num2cell(LinesEq)]; xlswrite(XLSInpFileName, LinesEqTxt, 'LinesEq')
end
