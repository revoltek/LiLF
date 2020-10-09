#!/bin/bash

# Exit on any error. More complex stuff could be done in future
# (see https://stackoverflow.com/questions/4381618/exit-a-script-on-error)
set -e

source /opt/lofar/lofarinit.sh

if [ "x$SAFE_MODE" == "xTrue" ]; then

    echo ""
    echo "[INFO] Not executing entrypoint as we are in safe mode, just opening a Bash shell."
    exec /bin/bash

else

    echo ""
    echo "[INFO] Executing entrypoint..."
    
    if [ "x$GUI" == "xTrue" ]; then
	    if [ "x$BASE_PORT" == "x" ]; then
	        echo "[INFO] No task base port set, will set noVNC port 8590 and VNC port 5900 with desktop id \"0\""  
	    else 
	        echo "[INFO] Task base port set, will set noVNC port $BASE_PORT and noVNC port $(($BASE_PORT+1)) with desktop id \"$(($BASE_PORT-5900+1))\""
	    fi
    fi
    
    #---------------------
    #   Setup home
    #---------------------

    if [ -f "/home/lofar/.initialized" ]; then
        :
    else
        echo "[INFO] Setting up home"
        mkdir -p /home/lofar

        # Copy over vanilla home contents
        for x in /lofar_home_vanilla/* /lofar_home_vanilla/.[!.]* /lofar_home_vanilla/..?*; do
            if [ -e "$x" ]; then cp -a "$x" /home/lofar/; fi
        done
        
        # Mark as initialized
        touch /home/lofar/.initialized
    fi
    

    #---------------------
    #   Save env
    #---------------------
    echo "[INFO] Dumping env"
    
    # Save env vars for later usage (e.g. ssh)
    
    env | \
    while read env_var; do
      if [[ $env_var == HOME\=* ]]; then
          : # Skip HOME var
      elif [[ $env_var == PWD\=* ]]; then
          : # Skip PWD var
      else
          echo "export $env_var" >> /tmp/env.sh
      fi
    done
    
    #---------------------
    #   VNC Password
    #---------------------
    if [ "x$GUI" == "xTrue" ]; then
	    if [ "x$AUTH_PASS" != "x" ]; then
	        echo "[INFO] Setting up VNC password..."
	        mkdir -p /home/lofar/.vnc
	        /opt/tigervnc/usr/bin/vncpasswd -f <<< $AUTH_PASS > /home/lofar/.vnc/passwd
	        chmod 600 /home/lofar/.vnc/passwd
	        export VNC_AUTH=True
	    else
	        echo "[INFO] Not setting up any VNC password"
	            
	    fi
    fi
    
	echo "[INFO] Creating /tmp/lofarhome to be used as lofar home"
	mkdir /tmp/lofarhome
	
	echo "[INFO] Initializing /tmp/lofarhome with configuration files"
	cp -aT /lofar_home_vanilla /tmp/lofarhome
	
	echo "[INFO] Moving to /home/lofar and setting as home"
	cd /home/lofar
	export HOME=/home/lofar
	
	echo "[INFO] Setting new prompt @$CONTAINER_NAME container"
	echo 'export PS1="${debian_chroot:+($debian_chroot)}\u@$CONTAINER_NAME@\h:\w\$ "' >> /tmp/lofarhome/.bashrc

	
    # Set entrypoint command
	if [ "x$@" == "x" ]; then
	    if [ "x$GUI" == "xTrue" ]; then
            COMMAND="supervisord -c /etc/supervisor/supervisord.conf"
	    else
	        COMMAND="/bin/bash"
	    fi
	else
	    COMMAND="$@"
	fi
	

    # Start!
	echo -n "[INFO] Will execute entrypoint command: "
	echo $COMMAND
	echo ""
	echo "=============================================================="
	echo ""
	echo "      Welcome to the LOFAR-IT $CONTAINER_NAME container!"
	echo ""
	echo "=============================================================="
	echo ""
	echo "You are now in /home/lofar with write access as user \"$(whoami)\"."
	echo ""
	echo "Remember that contents inside this container, unless stored"
	echo "on a persistent volume mounted from you host machine, will"
	echo "be wiped out when exiting the container."
	echo ""
	
	exec $COMMAND

fi
