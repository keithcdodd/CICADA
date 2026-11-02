function digest = cicada_task_glm_sha256(filePath)
% CICADA_TASK_GLM_SHA256  Stream a local input/source file as raw bytes.
% Read only; portable across MATLAB installations with Java enabled.
% Uses chunks so no full 4D BOLD is loaded by the hashing function.
assert(ischar(filePath) || (isstring(filePath) && isscalar(filePath)), ...
    'SHA256 path must be text');
filePath=char(filePath);
assert(isfile(filePath),'SHA256 input does not exist: %s',filePath);
fid=fopen(filePath,'rb');
assert(fid>=0,'Cannot open SHA256 input: %s',filePath);
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
% Read only the expected number of bytes. An extra fread past EOF can
% set MATLAB's EOF indicator and make a normal read appear to have failed.
fileInfo=dir(filePath);
assert(isscalar(fileInfo),'Cannot inspect SHA256 input size: %s',filePath);
remaining=double(fileInfo.bytes);
while remaining>0
    nToRead=min(65536,remaining);
    [chunk,nRead]=fread(fid,nToRead,'*uint8');
    assert(nRead==nToRead,'Short read while hashing %s',filePath);
    md.update(typecast(chunk,'int8'));
    remaining=remaining-nRead;
end
[msg,errnum]=ferror(fid);
assert(errnum==0,'File read failed while hashing %s: %s',filePath,msg);
raw=typecast(int8(md.digest()),'uint8');
digest=lower(reshape(dec2hex(raw,2).',1,[]));
assert(numel(digest)==64,'Unexpected SHA256 digest length');
end
