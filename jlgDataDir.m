function thePath = jlgDataDir

if exist('/jukebox','file')
    thePath = '/jukebox/tank/jlgauthi/Data';
else
    thePath = '/Volumes/tank/jlgauthi/Data';
end

