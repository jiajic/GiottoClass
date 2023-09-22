Giotto class diagram
================

## Class diagram and inheritance

``` mermaid
classDiagram

class nameData {
    [VIRTUAL]
    +name: character
}

class exprData {
    [VIRTUAL]
    +exprMat: ANY
}

class coordDataDT {
    [VIRTUAL]
    +coordinates: data.table
}

class metaData {
    [VIRTUAL]
    +metaDT: data.table
    +col_desc: character
}

class enrData {
    [VIRTUAL]
    +method: character
    +enrichDT: data.table
}

class nnData {
    [VIRTUAL]
    +nn_type: character
    +igraph: ANY
}

class spatNetData {
    [VIRTUAL]
    +method: character
    +parameters: ANY
    +outputObj: ANY
    +networkDT: data.table
    +networkDT_before_filter: data.table
    +cellShapeObj: ANY
}

class spatGridData {
    [VIRTUAL]
    +method: character
    +parameters: ANY
    +gridDT: data.table
}

class provData {
    [VIRTUAL]
    +provenance: ANY
}

class spatData {
    [VIRTUAL]
    +spat_unit: character
}

class featData {
    [VIRTUAL]
    +feat_data: character
}

class miscData {
    [VIRTUAL]
    +misc: ANY
}

class spatFeatData {
    [VIRTUAL]
}

class dimObj {
    +reduction: character
    +reduction_method: character
    +coordinates: ANY
    +misc ANY
}

class spatialNetworkObj {
    +crossSectionObjects: ANY
}

class giottoPolygon {
    +spatVector: ANY
    +spatVectorCentroids: ANY
    +overlaps: ANY
    +unique_ID_cache: character
}

class giottoPoints {
    +spatVector: ANY
    +networks: ANY
    +unique_ID_cache: character
}




provData --> spatData

spatData --> spatFeatData
featData --> spatFeatData

nameData --> exprObj
exprData --> exprObj
spatFeatData --> exprObj
miscData --> exprObj
giottoSubobject --> exprObj

metaData --> cellMetaObj
spatFeatData --> cellMetaObj
giottoSubobject --> cellMetaObj

metaData --> featMetaObj
spatFeatData --> featMetaObj
giottoSubobject --> featMetaObj

nameData --> dimObj
spatFeatData --> dimObj
giottoSubobject --> dimObj

nameData --> nnNetObj
nnData --> nnNetObj
spatFeatData --> nnNetObj
miscData --> nnNetObj
giottoSubobject --> nnNetObj

nameData --> spatLocsObj
coordDataDT --> spatLocsObj
spatData --> spatLocsObj
miscData --> spatLocsObj
giottoSubobject --> spatLocsObj

nameData --> spatialNetworkObj
spatNetData --> spatialNetworkObj
spatData --> spatialNetworkObj
miscData --> spatialNetworkObj
giottoSubobject --> spatialNetworkObj

nameData --> spatialGridObj
spatGridData --> spatialGridObj
spatFeatData --> spatialGridObj
miscData --> spatialGridObj
giottoSubobject --> spatialGridObj

nameData --> spatEnrObj
enrData --> spatEnrObj
spatFeatData --> spatEnrData
miscData --> spatEnrData
giottoSubobject --> spatEnrData

nameData --> giottoPolygon
giottoSubobject --> giottoPolygon

featData --> giottoPoints
giottoSubobject --> giottoPoints



```