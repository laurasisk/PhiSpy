import types
from copy import copy


class SeqioFilter( list ):
    """This class is to allow filtering of the Biopython SeqIO record

    SeqIO.parse returns a generator object so anytime you want to perform
    an action on it, you must iterate through the entire list. This
    class adds the ability to filter and return only a subset of the
    features. Those are split into separate features or joined based
    on the distance between the parts and then sorted based on the
    start to end location.

    Note:
        To use simply pass a SeqIO.parse object to it and then when
        the object is called a keyword is passed to it and only those
        features matching the keyword are returned.
    Example:
        record = SeqioFilter(SeqIO.parse(infile, infile_type))):
        #no change to standard SeqIO calls
        for entry in record:
            print(entry.id, entry.seq)
        #now we can get only certain features
        for cds in record.get_features('CDS'):
            print(cds)

    """

    def __init__( self, content ):
        self.n = 0
        self.get_n = dict()
        for n, item in enumerate(content):
            #for feature in item.features:
            #    if feature.location_operator == 'join':
            #        self.merge_or_split(feature)
            #item.features.sort( key = lambda feature : tuple([min(feature.location.start, feature.location.end),  max(feature.location.start, feature.location.end)]) )

            self.append(item)
            self.get_n[item.id] = n
            self.attach_methods(item)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            item = self[self.n]
        except IndexError:
            self.n = 0
            raise StopIteration()
        self.n += 1
        return item

    def __call__( self, keyword='' ):
        pass

    def get_entry(self, id):
        return self[self.get_n[id]]

    def distance_between(self, locations):
        return 111

    def merge_or_split(self, feature):
        cutoff_distance = 100
        feature.location_operator = None

        if self.distance_between(feature.location.parts) < cutoff_distance:
            #do merging stuff here
            return None
        else:
            #do splitting stuff here
            feature_copy = copy(feature)
            return feature_copy
    
 
    def attach_methods(self, target):
        """This method allows attaching new methods to the SeqIO entry object

           Args:
               target: is the SeqIO object that will be attaching a method to
        """

        # if feature has a complex location then merge or split it
        new_features = []
        for feature in target.features:
            if feature.location_operator == 'join':
                new_features.append(self.merge_or_split(feature))
        target.features.extend(new_features)

        # sort the features based on their location in the genome 
        target.features.sort( key = lambda feature : tuple([min(feature.location.start, feature.location.end),  max(feature.location.start, feature.location.end)]) )

        def get_features(target, feature_type):
            for feature in target.features:
                feature.id       = " ".join(feature.qualifiers.get('locus_tag', [str(feature.location)]))
                feature.function = " ".join(feature.qualifiers.get('product',['unknown']))
                # set start and stop for the feature and swap for direction if needed
                feature.start    = int(feature.location.start) + 1
                feature.stop     = int(feature.location.end)
                if feature.strand < 0:
                    feature.start, feature.stop = feature.stop, feature.start
                # if a feature type was provided, only return those features
                if not feature_type or feature.type == feature_type:
                    yield feature
        target.get_features = types.MethodType(get_features,target)

