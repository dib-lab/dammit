import pandas as pd

sources = ['HMMER',
           'LAST',
           'transdecoder',
           'crb-blast',
           'Infernal']

def create_tracks(transcript, transcript_df, track_start=100, track_width=40, track_sep=10):
    tracks = []
    for source in sources:
        if source in transcript_df.source.values:
            source_df = transcript_df.query('source == "{0}"'.format(source))
            track = {}
            track['trackName'] = source
            track['trackType'] = "track"
            track['trackFeatures'] = "complex"
            track['visible'] = True
            track['min_slice'] = True
            track['showTooltip'] = True

            track['inner_radius'] = track_start
            track['outer_radius'] = track_start + track_width

            items = []
            id_count = 0
            for _, row in source_df.iterrows():

                if row.notnull().Note:
                    name = row.Note
                elif row.notnull().Name:
                    name = row.Name
                else:
                    name = row.ID

                items.append({'id': id_count,
                              'start': int(row.start),
                              'end': int(row.end),
                              'name': name
                              })
                id_count += 1
            track['items'] = items
            tracks.append(track)
        track_start += track_width + track_sep

    return tracks
