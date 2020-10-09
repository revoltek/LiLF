# This is the Stager API wrapper module for the Lofar LTA staging service.
#
# It uses an xmlrpc proxy to talk and authenticate to the remote service. Your account credentials will be read from
# the awlofar catalog Environment.cfg, if present or can be provided in a .stagingrc file in your home directory.
#
# !! Please do not talk directly to the xmlrpc interface, but use this module to access the provided functionality.
# !! This is to ensure that when we change the remote interface, your scripts don't break and you will only have to
# !! upgrade this module.

__version__ = "1.3"

import datetime
from os.path import expanduser
from astropy.utils.compat.futures import process

# Python2/3 dependent stuff
from sys import version_info
python_version = version_info.major
if python_version == 3:
    import xmlrpc.client as xmlrpclib
    string_types = str
else:
    import xmlrpclib
    string_types = basestring

#---
# Determine credentials and create proxy
user = None
passw = None
try:
    f = expanduser("~/.awe/Environment.cfg")
    with open(f,'r') as file:
        print("%s - stager_access: Parsing user credentials from \"%s\"" % (datetime.datetime.now(), f))
        for line in file:
            if line.startswith("database_user"):
                user = line.split(':')[1].strip()
            if line.startswith("database_password"):
                passw = line.split(':')[1].strip()
except IOError:
    f = expanduser("~/.stagingrc")
    with open(f,'r') as file:
        print("%s - stager_access: Parsing user credentials from \"%s\"" % (datetime.datetime.now(), f))
        for line in file:
            if line.startswith("user"):
                user = line.split('=')[1].strip()
            if line.startswith("password"):
                passw = line.split('=', 1)[1].strip()

print("%s - stager_access: Creating proxy" % (datetime.datetime.now()))
proxy = xmlrpclib.ServerProxy("https://"+user+':'+passw+"@webportal.astron.nl/service-public/xmlrpc")
LtaStager =  proxy.LtaStager

# ---

def stage(surls):
    """ Stage list of SURLs """
    if isinstance(surls, str):
        surls = [surls]
    stageid = proxy.LtaStager .add_getid(surls)
    return stageid

def get_status(stageid):
    """ Get status of request with given ID """
    return proxy.LtaStager.getstatus(stageid)

def abort(stageid):
    """ Abort running request / release data of a finished request with given ID """
    return proxy.LtaStager.abort(stageid)

def get_surls_online(stageid):
    """ Get a list of all files that are already online for a running request with given ID  """
    return proxy.LtaStager.getstagedurls(stageid)

def get_srm_token(stageid):
    """ Get the SRM request token for direct interaction with the SRM site via Grid/SRM tools """
    return proxy.LtaStager.gettoken(stageid)

def reschedule(stageid):
    """ Reschedule a request with a given ID, e.g. after it was put on hold due to maintenance """
    return proxy.LtaStager.reschedule(stageid)

def get_progress(status=None, exclude=False):
    """ Get a detailed list of all running requests and their current progress.
        As a normal user, this only returns your own requests.
        :param status: If set to a valid status then only requests with that
        status are returned.
        :param exclude: If set to True then the requests with status 'status' are
        excluded.
    """
    all_requests = proxy.LtaStager.getprogress()
    if status is not None and isinstance(status, string_types):
        if python_version == 3:
            all_items = all_requests.items()
        else:
            all_items = all_requests.iteritems()
        if exclude is False:
            requests = {key: value for key, value in all_items if value["Status"] == status}
        else:
            requests = {key: value for key, value in all_items if value["Status"] != status}
    else:
        requests = all_requests
    return requests

def reschedule_on_status(status=None):
    """ Reschedule requests that have a status "on hold" or "aborted".
        :param status: The status that a request has to have in order to be
        rescheduled.
    """
    if status is not None and isinstance(status, string_types) and (status == "on hold" or status == "aborted"):
        requests = get_progress(status)
        for key in requests.keys():
            reschedule(int(key))
    else:
        print("The parameter status is either None, not a string neither of \"on hold\" nor \"aborted\".")

def get_storage_info():
    """ Get storage information of the different LTA sites, e.g. to check available disk pool space. Requires support role permissions. """
    return proxy.LtaStager.getsrmstorageinfo()

def prettyprint(dictionary, indent=""):
    """ Prints nested dict responses nicely. Example: 'stager_access.prettyprint(stager_access.get_progress())'"""
    if type(dictionary) is dict:
        for key in sorted(dictionary.keys()):
           item = dictionary.get(key)
           if type(item) is dict:
              print("%s+ %s" % (indent, str(key)))
              prettyprint(item, indent=indent+'  ')
           else:
              print("%s- %s\t->\t%s" % (indent, str(key), str(item)))
    else:
	    print("stager_access: This prettyprint takes a dict only!")

def reschedule_on_hold():
    """ Reschedule requests that are on hold.
    """
    reschedule_on_status("on hold")

def print_on_hold():
    """Print a list of all requests that are on hold.
    """
    requests = get_progress("on hold")
    prettyprint(requests)

def reschedule_aborted():
    """ Reschedule requests that were aborted.
    """
    reschedule_on_status("aborted")

def print_aborted():
    """Print a list of all requests that were aborted.
    """
    requests = get_progress("aborted")
    prettyprint(requests)

def print_running():
    """Print a list of requests that are currently executed and do not have a
    status of "success"
    """
    requests = get_progress("success", True)
    prettyprint(requests)
