import labkey

# classes to retrieve and add data to dholk.primate.wisc.edu LabKey Server

class LabkeyConnect:    # define Labkey server and project
    def serverContext(self, labkey_path):
        '''define server context by specifying server URL and path to folder to search'''
        self.labkey_server_url = 'dholk.primate.wisc.edu'
        self.labkey_folder_path  = labkey_path
        self.server_context = labkey.utils.create_server_context(self.labkey_server_url, self.labkey_folder_path, use_ssl=True, api_key='apikey|XXX')
        # self.server_context = labkey.utils.create_server_context(self.labkey_server_url, api_key='apikey|XXX', self.labkey_folder_path, use_ssl=True)
        # api_key='apikey|XXX',

class LabkeySelectRows(LabkeyConnect):

    # subclass for selecting rows that match a specific filter criteria
    # right now, allows only one filter criteria, though this could be expanded if necessary

    def set_filters(self, column, value):
        '''construct filter array for selecting rows'''
        self.labkey_filters = [labkey.query.QueryFilter(column, value)]

    def selectRows(self, labkey_schema, labkey_table):
        '''select data from LabKey Server'''
        self.labkey_schema =  labkey_schema
        self.labkey_table = labkey_table

        # query labkey
        result = labkey.query.select_rows(self.server_context,
                                      self.labkey_schema,
                                      self.labkey_table,
                                      filter_array=self.labkey_filters)

        return result

class LabkeyInsertRows(LabkeyConnect):

    # subclass for inserting rows to a LabKey table
    # right now, allows only one filter criteria, though this could be expanded if necessary

    def insertRows(self, labkey_schema, labkey_table, rows):
        '''select data from LabKey Server'''
        self.labkey_schema =  labkey_schema
        self.labkey_table = labkey_table

        # query labkey
        result = labkey.query.insert_rows(self.server_context,
                                      self.labkey_schema,
                                      self.labkey_table,
                                      rows)

        return result
