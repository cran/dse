\name{Sys.mail}
\alias{Sys.mail}
\title{Send Mail}
\usage{
Sys.mail(address= Sys.info()$user, ccaddress = NULL, bccaddress = NULL,
           subject = NULL, body= "no body", method = getOption("mailer"))
}
\arguments{
  \item{address}{Recipient's email address. The default is the user name as
     indicated by \code{Sys.info()$user}, which will correspond to your 
     own email address on many systems.}
  \item{ccaddress}{Optional email address for copies. 
    Use \code{ccaddress = NULL} for no copies.}
  \item{bccaddress}{Optional email address for blind copies.
    Use \code{bccaddress = NULL} for no copies.}
  \item{subject}{Subject of the email. Please do not use single quotes
     in the subject!}
  \item{method}{Submission method, one of \code{"mailx"},
    \code{"mail"}, or \code{"Mail"}.}
  \item{body}{A vector of character, one line per element to place in the body
     of the email.}
}
\description{
  Uses a system mail tool to send a message. Multiple recipients can be
  specified in any of the addresses by enclosing them in double (\code{"}) 
  quotes. Please do not use single (\code{'}) quotes.
  
  The body is a vector of character. Each element of body is written 
  as a line of the message body.
}
\details{
  Currently works only on Unix systems using "mailx", "mail" or "Mail". The mail
  tool needs to be found on your Unix search path.
}
\value{TRUE or FALSE indicating success or failure as returned by the operating
    system. (This does not indicate that the message was received at the 
    address.)}

\examples{
    \dontrun{Sys.mail(address="friends", subject="R", body="R is great")}
}
\seealso{\code{link{options}} }
\author{Paul Gilbert}
\keyword{utilities}


