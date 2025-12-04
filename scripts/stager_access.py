#! /usr/bin/env python3

# Copyright 2020 Stichting Nederlandse Wetenschappelijk Onderzoek Instituten,
# ASTRON Netherlands Institute for Radio Astronomy
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This is the Stager API wrapper module for the Lofar LTA staging service.
#
# It uses the stageit (lofar stager) api to talk to the remote service.
# Your api_token (see README) will be read from the awlofar catalog Environment.cfg, if present
# or can be provided in a .stagingrc file in your home directory.


__version__ = "2.0"

import datetime
import stager_proxy
proxy = stager_proxy.StagerProxy()

def deprecated(func):
    def inner():
        print(func.__name__, 'is deprecated. '
                             'Please consider reading the stageit documentation for usage')
        raise NotImplementedError('Use of deprecated functionality')

    return inner

def stage(surls, send_notifications=True):
    """ Stage list of SURLs, optionally enable/disable email notifications """
    if isinstance(surls, str):
        surls = [surls]
    try:
        return proxy.LtaStager.add_and_get_id(surls, send_notifications)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def get_status(stageid):
    """ Get status of request with given ID """
    try:
        return proxy.LtaStager.get_status(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def abort(stageid):
    """ Abort running request / release data of a finished request with given ID """
    try:
        return proxy.LtaStager.abort(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def get_surls_requested(stageid):
    """ Get a list of all files that are requested via a running request with given ID  """
    try:
        return proxy.LtaStager.get_requested_urls(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

# filter on response props (add to graphql)
def get_surls_pending(stageid):
    """ Get a list of all files that are not yet online for a running request with given ID  """
    try:
        return proxy.LtaStager.get_outstanding_urls(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

# filter on response props (add to graphql)
def get_surls_failed(stageid):
    """ Get a list of all files that have failed in a running request with given ID  """
    try:
        return proxy.LtaStager.get_failed_urls(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

# filter on response props (add to graphql)
def get_surls_online(stageid):
    """ Get a list of all files that are already online for a running request with given ID  """
    try:
        return proxy.LtaStager.get_staged_urls(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def get_srm_token(stageid):
    """ Get the SRM request token for direct interaction with the SRM site via Grid/SRM tools """
    try:
        return proxy.LtaStager.get_token(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def reschedule(stageid):
    """ Reschedule a request with a given ID, e.g. after it was put on hold due to maintenance """
    try:
        return proxy.LtaStager.reschedule(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def get_progress(status=None, exclude=False, only_last24h=False):
    """ Get a detailed list of all requests and their current progress.
        As a normal user, this only returns your own requests.
        :param status: If set to a valid status then only requests with that
        status are returned.
        :param exclude: If set to True then the requests with status 'status' are
        excluded.
        :param only_last24h: When True only requests that have been active in the last 24h get returned.  This can avoid the nasty time-outs that occasionally happen when the stager is stuck with too many retried requests.
    """
    since = datetime.datetime.now() - datetime.timedelta(hours=24) if only_last24h else datetime.datetime.min
    try:
        return proxy.LtaStager.get_progress(status=status, exclude=exclude, since=since)
    except Exception as e:
        proxy.print_exception_and_exit(e)

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

def print_aborted():
    """Print a list of all requests that were aborted.
    """
    requests = get_progress("A")
    prettyprint(requests)

def print_running():
    """Print a list of requests that are currently executed and do not have a
    status of "completed"
    """
    requests = get_progress("C", True)
    prettyprint(requests)

def get_macaroons(stageid):
    """ Get the macaroon(s) of the requested stageid. To be used when downloading data via webdav """
    try:
        return proxy.LtaStager.get_macaroons(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def get_webdav_urls_requested(stageid):
    """ Get the webdav urls of the requested stageid """
    try:
        return proxy.LtaStager.get_webdav_urls(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

def get_gridftp_urls_requested(stageid):
    """ Get the gridftp urls of the requested stageid"""
    try:
        return proxy.LtaStager.get_gridftp_urls(stageid)
    except Exception as e:
        proxy.print_exception_and_exit(e)

@deprecated
def get_stuck_requests(min_retries = 3, days_without_activity = 30):
    pass

@deprecated
def reschedule_on_status(status=None):
    pass

@deprecated
def get_storage_info():
    pass

@deprecated
def reschedule_on_hold():
    pass

@deprecated
def print_on_hold():
    pass

@deprecated
def reschedule_aborted():
    pass
