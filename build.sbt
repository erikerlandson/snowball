
// do clean first!
// sbt clean
// sbt unidoc
// sbt previewSite
// sbt ghpagesPushSite

// initiate publish from clean repo
// sbt clean
// sbt publish

name := "snowball"

version := "0.3.0"

organization := "com.manyangled"

//isSnapshot := true,

//publishConfiguration := publishConfiguration.value.withOverwrite(true)

//publishLocalConfiguration := publishLocalConfiguration.value.withOverwrite(true)

pomIncludeRepository := { _ => false }

publishMavenStyle := true

publishTo := {
  val nexus = "https://oss.sonatype.org/"
  if (isSnapshot.value)
    Some("snapshots" at nexus + "content/repositories/snapshots")
  else
    Some("releases"  at nexus + "service/local/staging/deploy/maven2")
}

licenses += ("Apache-2.0", url("http://opensource.org/licenses/Apache-2.0"))

homepage := Some(url("https://github.com/erikerlandson/snowball/"))

scmInfo := Some(
  ScmInfo(
    url("https://github.com/erikerlandson/snowball.git"),
    "scm:git@github.com:erikerlandson/snowball.git"
  )
)

developers := List(
  Developer(
    id    = "erikerlandson",
    name  = "Erik Erlandson",
    email = "eje@redhat.com",
    url   = url("https://erikerlandson.github.io/")
  )
)

crossPaths := false // drop off Scala suffix from artifact names.

autoScalaLibrary := false // exclude scala-library from dependencies

resolvers ++= Seq(
  Resolver.sonatypeRepo("releases"),
  Resolver.sonatypeRepo("snapshots")
)

// commons math used to be '% Provided' but the 'packageDoc' target
// now fails with that, so I'm just going to make it required
libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.6.1",
  "com.manyangled" % "gibbous" % "0.3.0",
  "com.manyangled" %% "gnuplot4s" % "0.2.0" % Test,
  "org.hamcrest" % "hamcrest-library" % "1.3" % Test,
  "com.novocode" % "junit-interface" % "0.11" % Test
)

compileOrder := CompileOrder.JavaThenScala

javacOptions ++= Seq()

enablePlugins(JavaUnidocPlugin, PublishJavadocPlugin, GhpagesPlugin)

siteSubdirName in JavaUnidoc := "java/api"

addMappingsToSiteDir(mappings in (JavaUnidoc, packageDoc), siteSubdirName in JavaUnidoc)

git.remoteRepo := "git@github.com:erikerlandson/snowball.git"
