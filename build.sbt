name := "snowball"

organization := "com.manyangled"

version := "0.1.0"

scalaVersion := "2.11.8"

crossScalaVersions := Seq("2.11.8", "2.12.4")

resolvers ++= Seq(
  Resolver.sonatypeRepo("releases"),
  Resolver.sonatypeRepo("snapshots")
)

libraryDependencies ++= Seq(
  //"com.joptimizer" % "joptimizer" % "4.0.0",
  "com.manyangled" %% "gibbous" % "0.1.0-local-d874ce55",
  "com.manyangled" %% "gnuplot4s" % "0.1.0-local-deedf561" % Test,
  "org.scalatest" %% "scalatest" % "3.0.5" % Test
)

licenses += ("Apache-2.0", url("http://opensource.org/licenses/Apache-2.0"))

scalacOptions ++= Seq("-unchecked", "-deprecation", "-feature")

scalacOptions in (Compile, doc) ++= Seq("-doc-root-content", baseDirectory.value+"/root-doc.txt")

enablePlugins(ScalaUnidocPlugin, GhpagesPlugin)

git.remoteRepo := "git@github.com:erikerlandson/snowball.git"
