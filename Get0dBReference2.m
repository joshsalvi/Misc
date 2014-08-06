% Returns the 0 dB reference value of the given signal
%
function db0 = Get0dBReference(file, signal)
%
if file.UseGlobalDbReferences == 1
    ptcPrefDbReferences = 0;
    dbRefs = get(file.preferences, 'Type', ptcPrefDbReferences);
else
    ptcInfoDbReferences = 11;
    dbRefs = get(file.infos, 'Type', ptcInfoDbReferences);
end
polymatlab = actxserver('PolyFile.PolyMatlab');
signal3 = invoke(polymatlab, 'Cast', 'ISignal3', signal);
db0 = invoke(signal3, 'Get0dB', dbRefs);
delete(polymatlab);
