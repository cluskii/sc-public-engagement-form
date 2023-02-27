from __future__ import print_function

import string

from googleapiclient.discovery import build


def extract_from_google_sheets(creds, sheet_id):
    service = build('sheets', 'v4', credentials=creds)

    # Call the Sheets API
    sheet = service.spreadsheets()
    result = sheet.values().get(spreadsheetId=sheet_id,
                                range="A:ZZ").execute()
    values = result.get('values', [])

    if not values:
        print('No data found.')
    else:
        return values

