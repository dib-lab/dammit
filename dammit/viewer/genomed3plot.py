import pandas as pd


def create_tracks(transcript, df, sources, track_start=100, track_width=40, track_sep=10,
                  trackname_func=lambda name: name.partition('.')[2]):

    transcript_df = df.query('seqid == "{0}"'.format(transcript))
    tracks = []
    for src, source_df in [(src,df.query('source == "{0}"'.format(src))) for src in sources]:
        if src in transcript_df.source.values:
            track = {}
            track['trackName'] = trackname_func(src) if trackname_func is not None else src
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
