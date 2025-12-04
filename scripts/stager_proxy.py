import datetime
import json
from os.path import expanduser, exists
from sys import version_info

from requests import HTTPError

try:
    import requests
except ImportError:
    print("ERROR: package 'requests' not found. Make sure you have the 'requests' package installed (or run 'pip install -r requirements.txt'")

python_version = version_info.major
if python_version == 3:
    import configparser
else:
    import ConfigParser as configparser
    import StringIO

STAGEIT_API = "stageit/api/staging"


def read_config_string(section, config_string):
    config = configparser.ConfigParser()
    section_and_config = '[{0}]\n{1}'.format(section, config_string)

    if python_version == 3:
        config.read_string(section_and_config)
        return config[section]
    else:
        buf = StringIO.StringIO(section_and_config)
        config.readfp(buf)
        return dict(config.items(section))


class StagerProxy:

    def __init__(self):
        self.proxy = None

    @staticmethod
    def print_exception_and_exit(e):
        if isinstance(e, HTTPError):
            print("Invalid HTTP response, something is wrong with the server")
        elif isinstance(e, requests.exceptions.Timeout):
            print("Timeout occurred, please try again")
        elif isinstance(e, requests.exceptions.TooManyRedirects):
            print("Too many redirects, did you change the api URL?")
        elif isinstance(e, requests.exceptions.ConnectionError):
            print("There is a network problem")
        elif isinstance(e, requests.exceptions.RequestException):
            print("General request error occurred")
        else:
            print("Unknown error: %s", str(e))

        raise SystemExit(e)


    @property
    def LtaStager(self):
        if self.proxy is None:
            self.proxy = self.create_proxy()
        return self.proxy

    def create_proxy(self):
        user, password, api_token, hostname = self.parse_credentials()
        # For test purposes you can add the next line in the .stagingrc file, it will then not use the production StageIT interface
        # hostname=sdc-dev.astron.nl
        if hostname:
            print(f"Use alternative hostname={hostname}")
        else:
            hostname = "sdc.astron.nl"
        if user is not None or password is not None:
            error = 'Found old stager user credentials, please acquire an access token from the new stager api (see README) and remove these old credentials'
            print(error)
            return None
        elif api_token is not None:
            print("%s - stager_access: Found API token. Creating json api proxy" % (datetime.datetime.now()))
            proxy = self.from_json_api(graphql_uri=f'https://{hostname}/{STAGEIT_API}/graphql/',
                                       rest_uri=f'https://{hostname}/{STAGEIT_API}/requests/',
                                       token=api_token)
            return proxy
        else:
            print("%s - stager_access:  Could not find a file with user credentials." % (datetime.datetime.now()))

    def parse_credentials(self):
        # Iterate over possible files that could contain user credentials.
        for f in ["~/.stagingrc"]:
            file = expanduser(f)
            try:
                credentials = self.parse_config_file(file)
                return credentials.get("user"), credentials.get("password"), credentials.get("api_token"), credentials.get("hostname"),
            except Exception:
                print("%s - stager_access: Could not parse user credential file '%s'." % (datetime.datetime.now(), file))

    def parse_config_file(self, file):
        """ Parse an ini file that shall contain user credentials and/or an api token
            for the LOFAR stager API.
        """
        if exists(file) is False:
            # File not found, nothing to do.
            raise Exception
        with open(file, 'r') as stream:
            print(
                "%s - stager_access: Parsing user credentials found in file '%s'." % (datetime.datetime.now(), file))

            # This may seem odd but the Python configparser cannot parse ini
            # files that just contain "A = B" lines without a section name
            # like "[Foo Bar]".  So I just add a section name just before the
            # config file is read.
            return read_config_string("LOFAR stager credentials", stream.read())

    @staticmethod
    def from_json_api(graphql_uri, rest_uri, token):
        return JsonApiStager(graphql_uri, rest_uri, token)


class JsonApiStager:
    def __init__(self, graphql_uri, rest_uri, token):
        self.graphql_uri = graphql_uri
        self.rest_uri = rest_uri
        self.token = token
        self.headers = {
            'Authorization': 'Token {0}'.format(token),
            'Content-Type': 'application/json'
        }

        self.get_request_graphql_query = """
{{
  request(id: {0}) {{
    {1}
  }}
}}
"""
        self.get_request_with_surls_graphql_query = """
{{
  request(id: {0}) {{
    surls {{
        surl
    }}
  }}
}}        
"""
        self.get_request_with_macaroons_graphql_query = """
        {{
          request(id: {0}) {{
            macaroons {{
                content,
                validUntil,
                ltaSite {{
                    name
                }}
            }}
          }}
        }}        
        """

        self.get_requests_graphql_query = """
{{
  requests(when_Gte:"{0}"{1}) {{
    edges {{
    	node {{
            id,
            requestToken,
            requestType,
            when,
            response,
            currentStatusCode,
            surls {{
              surl
            }}
      }}
    }}
  }}
}}        
"""

    def graphql_query(self, query, *keys):
        response = requests.post('{0}'.format(self.graphql_uri), json={'query': query}, headers=self.headers)
        response.raise_for_status()
        json = response.json()
        errors = json.get('errors', {})
        if len(errors) > 0:
            print(errors[0].get('message', 'unknown error'))
        else:
            return self.get_from_json(json, *keys)

    def rest_api_call(self, method, endpoint, data):
        response = requests.request(method, '{0}{1}'.format(self.rest_uri, endpoint), json=data, headers=self.headers)
        response.raise_for_status()

        # stageit backend returns 202 for empty bodies
        if response.status_code != 202:
            return response.json()

    @staticmethod
    def get_from_json(json, *keys):
        for key in keys:
            if isinstance(json, dict):
                json = json.get(key)
            else:
                return None

        return json

    def add_and_get_id(self, surls, send_notifications):
        json = self.rest_api_call('post', '', { 'surls': surls, 'notify': send_notifications, 'request_type': 'stage' })
        return json.get('id')

    def get_status(self, stageid):
        field = 'currentStatus'
        query = self.get_request_graphql_query.format(stageid, field)
        return self.graphql_query(query, 'data', 'request', field)

    def get_requested_urls(self, stageid):
        query = self.get_request_with_surls_graphql_query.format(stageid)
        surls = self.graphql_query(query, 'data', 'request', 'surls')
        return [surl.get('surl', '') for surl in surls]

    def get_outstanding_urls(self, stageid):
        lta_request_response = self._get_lta_request_response(stageid)
        return lta_request_response.get('progressing', [])

    def get_failed_urls(self, stageid):
        lta_request_response = self._get_lta_request_response(stageid)
        response_errors = lta_request_response.get('errors', {})
        return list(response_errors.keys())

    def get_staged_urls(self, stageid):
        lta_request_response = self._get_lta_request_response(stageid)
        return lta_request_response.get('online', [])

    def get_token(self, stageid):
        field = 'requestToken'
        query = self.get_request_graphql_query.format(stageid, field)
        return self.graphql_query(query, 'data', 'request', field)

    def abort(self, stageid):
        self.rest_api_call('put', str(stageid) + '/abort/', {})
        return True

    def reschedule(self, stageid):
        self.rest_api_call('patch', str(stageid) + '/', { 'current_status_code': 'S' })
        return 'scheduled'

    def get_progress(self, status=None, exclude=False, since=datetime.datetime.min):
        status_filter = ""

        if status is not None:
            status_filter = ", currentStatusCode_Ne: {0}" if exclude else ", currentStatusCode: {0}"
            status_filter = status_filter.format(status)

        query = self.get_requests_graphql_query.format(since.isoformat(), status_filter)
        nodes = self.graphql_query(query, 'data', 'requests', 'edges')
        requests = [node.get('node', {}) for node in nodes]
        requests = { r['id'] : r for r in requests }
        return requests

    def _get_lta_request_response(self, stageid):
        field = 'response'
        query = self.get_request_graphql_query.format(stageid, field)
        lta_request_response_json = self.graphql_query(query, 'data', 'request', field)
        lta_request_response = json.loads(lta_request_response_json)
        return lta_request_response

    def get_macaroons(self, stageid):
        query = self.get_request_with_macaroons_graphql_query.format(stageid)
        macaroons_result = self.graphql_query(query, 'data', 'request', 'macaroons')
        return [{ mac.get('ltaSite',).get('name', ''): mac.get('content', '')} for mac in macaroons_result]

    def get_webdav_urls(self, stageid):
        return self.rest_api_call('get', str(stageid) + '/convert2webdav/', {})

    def get_gridftp_urls(self, stageid):
        return self.rest_api_call('get', str(stageid) + '/convert2gridftp/', {})
