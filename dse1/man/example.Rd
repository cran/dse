\name{example}
\alias{example.BOC.93.4.paper.tests}
\alias{example.convert.eg.data}
\alias{example.get.eg.raw.data}
\alias{example.like.sub.sample}
\alias{example.like.total.sample}
\alias{example.show.ytoy.cpi}
\alias{example.tests}
\alias{example.truncate.data}
\alias{example.VAR.SVD}
\alias{example.verify.data}
\alias{example.verify.raw.data}

\title{TS Examples}
\description{DSE Package Examples}
\usage{
    example.BOC.93.4.paper.tests(example.data, example.raw.data,
       verbose=T, fuzz.small=1e-10)
    example.convert.eg.data(example.raw.data)
    example.get.eg.raw.data(file)
    example.like.sub.sample(model)
    example.like.total.sample(model)
    example.show.ytoy.cpi(model, raw.data, start = 1)
    example.tests(example.data, verbose=T, summary=T,
        fuzz.small=1e-14, fuzz.large=1e-8)
    example.truncate.data(d.all)
    example.VAR.SVD(d,d.all, d.all.raw)
    example.verify.data(example.data, verbose=T, fuzz.small=1e-14)
    example.verify.raw.data(example.raw.data, verbose=T, fuzz.small=1e-5)
}

\details{
    These functions run various examples. Their purpose is mainly to 
    illustrate the use of other functions. Print out the functions and
    examine them.
}

%\keyword{DSE}
\keyword{ts}

