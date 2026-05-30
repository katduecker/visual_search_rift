function [] = NAtA_bbox()

KbName('UnifyKeyNames'); % for easy use of Keyboard keys

activeKeysOnly=0; %If you want to listen only to specific keys, you can
device=0; %default device

if activeKeysOnly %here we listen only to specific keys
    active=[KbName('4$') KbName('7&')]; %These are the left and right index fingers of the (5-button) NATA boxes
    keylist=zeros(1,256); %Set all keys to zero (ignore)
    keylist(active)=1; %set active keys to 1 (listen)
    KbQueueCreate(device,keylist);%%Create queue, this is a time consuming operation (relatively), do while non-time critical
    disp('I will only listen to: ');disp(KbName(active));
else
    KbQueueCreate; %Default initialization, all keys listened to
end

noResponse=1; %While loop break

KbQueueStart(); %Start listening
startTime = GetSecs; %to get reaction times relative to some event (e.g. stimulus onset)

KbQueueFlush();% clear all keyboard presses so far. Basically, you want the responses after this point 

%listen to reponses for a fixed amount of time, using KbQueueCheck (similar to KbCheck)
disp('Press some keys for the next 10 seconds..');
while noResponse
    
    [pressed, firstpress]=KbQueueCheck(); %check response, return whether pressed, and first press timestamp
    
    %Note that two keys may have been pressed    
    keyCode=find(firstpress);

    if length(keyCode)>1 %two or more buttons pressed
        [~,ind]=min(firstpress(keyCode));
        keyCode = keyCode(ind); %select first response
    end

    t_keypress=firstpress(keyCode);
    
    if pressed
        %disp(['You pressed ' KbName(firstpress) ' ' num2str((firstpress(find(firstpress))-startTime)*1000) 'ms after start'])        
        disp(['You pressed ' KbName(keyCode) ' ' num2str((t_keypress-startTime)*1000) 'ms after start'])        
    end
    
    if (GetSecs-startTime)>10 %let's stop listening after 10 seconds
        noResponse=0;
    end
    
    WaitSecs(0.01); %no need to poll very tightly as responses are timestamped as they come in (unless you want something to happend right after the button press).
    
end

disp('Ok, I stopped listening.')

KbQueueFlush();% clear all keyboard presses so far. Basically, you want the responses after this point 

%now wait for a response (similar to KbWait)
disp('Press a key to continue...')
KbQueueWait();

disp('Great, done!')

KbQueueRelease() %clean up the KbQueue, you need to re-create in order to use it again.


end